/***
 *
 * Online Dynamic Connectivity (Holm, de Lichtenberg, Thorup)
 * Connectivity, component sizes and component value sums of an undirected graph under edge insertions and deletions
 *
 * Complexity:
 *   - add_edge, remove_edge: amortized O(log^2 n)
 *   - every other call: amortized O(log n)
 *   - memory: O((n + m) log n) for m add_edge calls, 420 MB for 1.3e6 calls churning a grid of 3e5 vertices
 *
 * OnlineDynamicConnectivity g(n): an empty graph on vertices 0 to n - 1, every value 0
 * OnlineDynamicConnectivity g(values): an empty graph on vertices 0 to n - 1 with the given values
 * g.add(v, x): adds x to the value of v
 * g.add_edge(u, v): adds a copy of the edge (u, v), parallel edges and self loops allowed
 * g.remove_edge(u, v): deletes one live copy of (u, v) in either orientation, the copy must exist
 * g.connected(u, v): u and v are in the same connected component
 * g.size(v): number of vertices in the component of v
 * g.sum(v): sum of the values in the component of v
 * g.components: number of connected components
 *
 * Notes:
 *   - Fully online; when all operations are known in advance, offline_dynamic_connectivity.cpp is simpler and faster
 *   - Everything is iterative, deep trees are safe
 *   - Algorithm from J. Holm, K. de Lichtenberg, M. Thorup, J. ACM 48(4), 2001, with Euler tour trees on splay trees
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct OnlineDynamicConnectivity{
    enum Mark : unsigned char{TREE = 1, NONTREE = 2};

    /// Level i keeps a spanning forest of the edges of level >= i, each tree as its cyclic Euler tour in one splay tree
    struct Node{
        long long sum;  /// first, so the node packs into 32 bytes
        int ch[2], par;
        int id;  /// a vertex slot if >= 0, else one of the two arcs of edge ~id
        int size;
        unsigned char mark, any;  /// TREE on one arc of each tree edge of this level, NONTREE on a vertex with a list
    };

    struct Edge{
        int u, v, level;  /// level -1 once deleted
        vector<int> arcs;  /// arcs[i] and arcs[i] + 1 are its arcs in level i's forest
    };

    /// Salted splitmix64: std::hash is the identity, so keys u * n + v chosen by a test can all share one bucket
    struct SplitMix{
        uint64_t salt;

        SplitMix() : salt(chrono::steady_clock::now().time_since_epoch().count()) {}

        /// noexcept, or unordered_map caches every hash code: 11-25% slower in euler_tour_tree
        size_t operator()(uint64_t x) const noexcept{
            x += salt + 0x9e3779b97f4a7c15;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
            x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
            return x ^ (x >> 31);
        }
    };

    int n, components;
    vector<long long> value;
    vector<Node> t;  /// nodes 0..n - 1 are the vertices of level 0, the only nodes holding values
    vector<int> vertex_of;  /// per vertex slot, its vertex; slot k is the k-th vertex node made, slot v < n is node v
    vector<vector<int>> nontree;  /// per vertex slot, non-tree edges at its level; stale once raised, linked or deleted
    vector<vector<int>> vertex_at;  /// per level, the node of each vertex, -1 until it gets a tree edge there
    vector<Edge> edges;
    vector<int> spare_arcs;
    unordered_map<long long, vector<int>, SplitMix> live;

    OnlineDynamicConnectivity(int n) : OnlineDynamicConnectivity(vector<long long>(n)) {}

    OnlineDynamicConnectivity(const vector<long long>& values)
        : n(values.size()), components(n), value(values), vertex_at(1, vector<int>(n, -1)){
        for (int v = 0; v < n; v++) node(0, v);
    }

    void add(int v, long long x){
        splay(v);
        value[v] += x;
        pull(v);
    }

    void add_edge(int u, int v){
        int e = edges.size();
        edges.push_back({u, v, 0, {}});
        live[key(u, v)].push_back(e);
        if (u == v) return;

        if (same_tree(u, v)) add_nontree(e);
        else{
            link(e);
            components--;
        }
    }

    void remove_edge(int u, int v){
        auto it = live.find(key(u, v));
        assert(it != live.end());
        int e = it->second.back();
        it->second.pop_back();
        if (it->second.empty()) live.erase(it);

        int level = edges[e].level;
        edges[e].level = -1;
        if (edges[e].arcs.empty()) return;

        for (int a: edges[e].arcs) cut(a);
        edges[e].arcs.clear();
        for (int i = level; i >= 0; i--){
            if (replace(i, u, v)) return;
        }
        components++;
    }

    bool connected(int u, int v){
        return same_tree(u, v);
    }

    int size(int v){
        splay(v);
        return t[v].size;
    }

    long long sum(int v){
        splay(v);
        return t[v].sum;
    }

    int node(int i, int v){
        if (i == (int)vertex_at.size()) vertex_at.emplace_back(n, -1);
        if (vertex_at[i][v] == -1){
            vertex_at[i][v] = t.size();
            t.push_back({i == 0 ? value[v] : 0, {-1, -1}, -1, (int)vertex_of.size(), 1, 0, 0});
            vertex_of.push_back(v);
            nontree.emplace_back();
        }

        return vertex_at[i][v];
    }

    long long key(int u, int v) const{
        return (long long)min(u, v) * n + max(u, v);
    }

    bool same_tree(int a, int b){
        splay(a);
        splay(b);
        return a == b || t[a].par != -1;
    }

    void add_nontree(int e){
        int i = edges[e].level;
        for (int w: {edges[e].u, edges[e].v}){
            int x = node(i, w);
            nontree[t[x].id].push_back(e);
            set_mark(x, NONTREE);
        }
    }

    /// Splices e into the forests of levels arcs.size()..level, its endpoints must be in different trees there
    void link(int e){
        while ((int)edges[e].arcs.size() <= edges[e].level){
            int i = edges[e].arcs.size(), a = new_arcs(e, i == edges[e].level ? TREE : 0);
            edges[e].arcs.push_back(a);

            /// v's tour A v B becomes A v (v->u) (u's tour from u) (u->v) B
            int ru = make_first(node(i, edges[e].u)), v = node(i, edges[e].v);
            splay(v);
            int after = t[v].ch[1];
            if (after != -1) t[after].par = -1;
            t[v].ch[1] = -1;
            join(join(join(v, a), ru), join(a + 1, after));
        }
    }

    void cut(int a){
        /// The two arcs split the cyclic tour into the part between them and two outer parts that close into one tour
        splay(a);
        splay(a + 1, a);
        int d = t[a].ch[1] == a + 1;
        int between = t[a + 1].ch[!d], outer_a = t[a].ch[!d], outer_b = t[a + 1].ch[d];
        for (int c: {between, outer_a, outer_b}){
            if (c != -1) t[c].par = -1;
        }
        join(outer_a, outer_b);
        spare_arcs.push_back(a);
    }

    /// After cutting a tree edge u-v of level >= i: searches level i for a replacement, raising the misses to i + 1
    bool replace(int i, int u, int v){
        int a = node(i, u), b = node(i, v);
        splay(a);
        splay(b);
        if (t[a].size > t[b].size) swap(a, b);

        bool raised = false;
        for (int x; (x = marked(a, NONTREE)) != -1; a = x){
            int slot = t[x].id, found = -1;
            while (found == -1 && !nontree[slot].empty()){
                int e = nontree[slot].back();
                nontree[slot].pop_back();
                if (edges[e].level != i || !edges[e].arcs.empty()) continue;

                int y = node(i, edges[e].u ^ edges[e].v ^ vertex_of[slot]);
                if (!same_tree(x, y)) found = e;
                else{
                    /// Non-tree edges of level i + 1 need their endpoints joined at level i + 1, so tree edges go first
                    if (!raised) raise_tree_edges(x);
                    raised = true;
                    edges[e].level++;
                    add_nontree(e);
                }
            }
            set_mark(x, nontree[slot].empty() ? 0 : NONTREE);

            if (found != -1){
                link(found);
                return true;
            }
        }

        return false;
    }

    /// Raising only the smaller tree's edges keeps every tree of level i + 1 at most half a tree of level i
    void raise_tree_edges(int a){
        for (int x; (x = marked(a, TREE)) != -1; a = x){
            set_mark(x, 0);
            int e = ~t[x].id;
            edges[e].level++;
            link(e);
        }
    }

    /// A node of x's tree with the mark bit, splayed to the root, or -1
    int marked(int x, int bit){
        splay(x);
        if (!(t[x].any & bit)) return -1;

        while (!(t[x].mark & bit)){
            int l = t[x].ch[0];
            x = l != -1 && (t[l].any & bit) ? l : t[x].ch[1];
        }
        splay(x);
        return x;
    }

    void set_mark(int x, unsigned char mark){
        splay(x);
        t[x].mark = mark;
        pull(x);
    }

    int new_arcs(int e, unsigned char mark){
        int a = t.size();
        if (!spare_arcs.empty()){
            a = spare_arcs.back();
            spare_arcs.pop_back();
        }
        else{
            t.resize(a + 2);
        }

        t[a] = {0, {-1, -1}, -1, ~e, 0, mark, mark};
        t[a + 1] = {0, {-1, -1}, -1, ~e, 0, 0, 0};
        return a;
    }

    /// The tour is cyclic, so rotating it to start at any node keeps it a valid tour; returns the new root
    int make_first(int x){
        splay(x);
        int left = t[x].ch[0];
        if (left == -1) return x;

        t[left].par = -1;
        t[x].ch[0] = -1;
        pull(x);
        return join(x, left);
    }

    /// Joins the tours rooted at a and b, a first; returns the new root
    int join(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        int x = a;
        while (t[x].ch[1] != -1) x = t[x].ch[1];

        splay(x);
        t[x].ch[1] = b;
        t[b].par = x;
        pull(x);
        return x;
    }

    void splay(int x, int goal = -1){
        while (t[x].par != goal){
            int p = t[x].par, g = t[p].par;
            if (g != goal) rotate((t[g].ch[1] == p) == (t[p].ch[1] == x) ? p : x);
            rotate(x);
        }
        pull(x);
    }

    /// Pulls only p, x stays stale until splay pulls it once at the end
    void rotate(int x){
        int p = t[x].par, g = t[p].par, d = t[p].ch[1] == x, child = t[x].ch[!d];

        if (g != -1) t[g].ch[t[g].ch[1] == p] = x;
        t[x].par = g;
        t[p].ch[d] = child;
        if (child != -1) t[child].par = p;
        t[x].ch[!d] = p;
        t[p].par = x;

        pull(p);
    }

    void pull(int x){
        Node& c = t[x];
        c.size = c.id >= 0;
        c.sum = x < n ? value[x] : 0;
        c.any = c.mark;
        for (int y: c.ch){
            if (y == -1) continue;
            c.size += t[y].size;
            c.sum += t[y].sum;
            c.any |= t[y].any;
        }
    }
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout, about 5x faster than scanf for integers, interactive problems included
 *
 * Complexity: O(1) amortized per character, one system call per buffer of up to 8192 bytes
 *
 * fio::flush(): writes the pending output now, which also happens at normal exit
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ..., returns the count read before the input ended
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - Reading into a char array does not compile: read a string, or a grid as vector<string>
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A - not followed by a digit reads as 0 and eats the byte after it, so "3 - 4" reads 3, then 0, then 4
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
 *   - A float, double or long double read gives strtof's, strtod's or strtold's value, 2-3x faster on short decimals
 *   - Interactive problems need no flush(): pending output is written each time the input buffer runs empty
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>
#include <unistd.h>  /// read() and write(), which MinGW's unistd.h also declares, so Windows judges compile them too

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// write, not fwrite: fwrite strands bytes in stdout's buffer, and fflush after it costs a second call, 7% slower
    void flush(){
        int done = 0;
        while (done < outptr){
            int len = ::write(1, outbuf + done, outptr - done);
            if (len <= 0) break;
            done += len;
        }
        outptr = 0;
    }

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            /// An interactive judge answers only after it has read our query, so the query leaves before we wait
            flush();

            /// read returns the bytes available, fread waits for all 8192 or EOF: a deadlock on an interactive pipe
            inptr = 0, buf_len = ::read(0, inbuf, BUF_SIZE);
            if (buf_len <= 0) return EOF;
        }

        return (unsigned char)inbuf[inptr++];
    }

    /// Catches every integer the overloads below miss, so __int128 needs no is_integral specialization (breaks C++20)
    template <typename T>
    bool read_one(T &x){
        int c = ' ', neg = 0;

        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// A local, since an int x may alias inptr or buf_len and writing it per digit forces their reload from memory
        T v = 0;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        while (isdigit(c)){
            v = v * 10 - (c - '0');

            /// The digits already buffered are read through a local pointer, 15% faster than a read_char call per digit
            const char *p = inbuf + inptr, *end = inbuf + buf_len;
            while (p < end && (unsigned)(*p - '0') < 10) v = v * 10 - (*p++ - '0');
            inptr = p - inbuf;

            /// A buffered non-digit ends the number and is eaten, else the digits may go on after a refill
            if (p < end){
                inptr++;
                break;
            }
            c = read_char();
        }
        x = neg ? v : -v;

        return true;
    }

    /// The C locale's isspace set, false for EOF, and 1.4-1.7x faster for strings and chars than isspace's table lookup
    inline bool is_space(int c){
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (is_space(c)) c = read_char();
        if (c == EOF) return 0;

        for (; !is_space(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return len;
    }

    /// Digits m and 10^|e| both exact in F round once, as parse does: 2-3x faster than parse on short tokens
    template <typename F, F (*parse)(const char *, char **)>
    bool read_floating(F &x, int max_len, unsigned long long max_mantissa, int max_exp){
        int len = read_token();
        if (!len) return false;

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds float and double twice, long double once
        bool own_width = !FLT_EVAL_METHOD || is_same<F, long double>::value;

        /// A longer token skips the digit scan, which would make its parse about 10% slower
        if (len <= max_len && own_width){
            static const F p10[] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L,
                                    1e10L, 1e11L, 1e12L, 1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L,
                                    1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};
            const unsigned char *s = (const unsigned char *)tmpbuf + (*tmpbuf == '-' || *tmpbuf == '+');
            unsigned long long m = 0;
            int digits = 0, e = 0;
            for (; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++;
            if (*s == '.'){
                for (s++; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++, e--;
            }

            /// No exponent digits after the e equals strtod, which then ignores the e
            if (*s == 'e' || *s == 'E'){
                int sign = s[1] == '-' ? -1 : 1, ex = 0;
                s += 1 + (s[1] == '-' || s[1] == '+');
                for (; isdigit(*s) && ex < 1000; s++) ex = ex * 10 + (*s - '0');
                e += sign * ex;
            }

            if (!*s && digits && digits <= 19 && m <= max_mantissa && e >= -max_exp && e <= max_exp){
                F v = e < 0 ? m / p10[-e] : m * p10[e];
                x = *tmpbuf == '-' ? -v : v;
                return true;
            }
        }

        x = parse(tmpbuf, nullptr);
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        return read_floating<float, strtof>(x, 10, 1 << 24, 10);
    }

    /// 22 characters keep 16-digit scientific tokens like -1.234567890123456e-05 on the fast path
    bool read_one(double &x){
        return read_floating<double, strtod>(x, 22, 1ULL << 53, 22);
    }

    /// Only an x87 64-bit mantissa takes the fast path: a double-double long double does not round a division once
    bool read_one(long double &x){
        return read_floating<long double, strtold>(x, 21, ~0ULL, LDBL_MANT_DIG == 64 ? 27 : -1);
    }

    bool read_one(string &s){
        int c = ' ';

        while (is_space(c)) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !is_space(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (is_space(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T> &v){
        for (auto &&x: v){
            if (!read_one(x)) return false;
        }

        return true;
    }

    template <typename T>
    int read(T &x){
        return read_one(x);
    }

    template <typename T, typename ...Args>
    int read(T &x, Args &...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
        return true;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(char c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    /// Integers, like read_one: the plain overloads below win for every other type
    template <typename T>
    void write_one(T x){
        int l = 0;

        if (x < 0){
            /// Peel one digit while still negative so that -x cannot overflow for the minimum value
            write_char('-');
            tmpbuf[l++] = '0' - x % 10;
            x = -(x / 10);
        }

        while (x || !l){
            tmpbuf[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmpbuf[--l]);
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char *s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    /// Without it a char * picks the integer template, an exact match that beats the const char * overload
    void write_one(char *s){
        write_one((const char *)s);
    }

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
        static const long long p10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
            10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 1000000000000000,
            10000000000000000, 100000000000000000, 1000000000000000000};

        if (f.digits >= 0 && f.digits <= 18){
            long double scaled = fabsl(f.value) * p10[f.digits];

            /// Rounding never moves the product across k + 0.5, exact below 2^52, so only an exact k + 0.5 needs printf
            if (scaled < 4503599627370496.0L && scaled - (long long)scaled != 0.5L){
                long long q = scaled + 0.5L;
                if (signbit(f.value)) write_char('-');
                write_one(q / p10[f.digits]);
                if (!f.digits) return;

                long long r = q % p10[f.digits];
                for (int j = f.digits - 1; j >= 0; j--) tmpbuf[j] = '0' + r % 10, r /= 10;
                write_char('.');
                for (int j = 0; j < f.digits; j++) write_char(tmpbuf[j]);
                return;
            }
        }

        int len = snprintf(tmpbuf, sizeof(tmpbuf), "%.*Lf", f.digits, f.value);

        if (len < (int)sizeof(tmpbuf)){
            for (int j = 0; j < len; j++) write_char(tmpbuf[j]);
            return;
        }

        string s(len + 1, '\0');
        snprintf(&s[0], len + 1, "%.*Lf", f.digits, f.value);
        s.resize(len);
        write_one(s);
    }

    template <typename T>
    void write_one(const vector<T> &v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    template <typename T>
    char separator_after(const T &){
        return ' ';
    }

    template <typename T>
    char separator_after(const vector<T> &){
        return '\n';
    }

    template <typename T>
    void write(const T &x){
        write_one(x);
        write_char('\n');
    }

    template <typename T, typename ...Args>
    void write(const T &x, const Args &...args){
        write_one(x);
        write_char(separator_after(x));
        write(args...);
    }

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }
}

using namespace fio;

int main(){
    int n, m, op, u, v, last = 0;

    read(n, m);
    OnlineDynamicConnectivity g(n + 1);
    while (m--){
        read(op, u, v);
        u ^= last, v ^= last;

        if (op == 0) g.add_edge(u, v);
        else if (op == 1) g.remove_edge(u, v);
        else{
            bool flag = g.connected(u, v);
            last = flag ? u : v;
            write(flag ? 'Y' : 'N');
        }
    }

    return 0;
}
