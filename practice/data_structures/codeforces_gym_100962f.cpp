/***
 *
 * Mo's Algorithm On Tree Paths
 * Offline path queries on a tree, by Mo's algorithm over an Euler tour that lists each node twice
 *
 * Complexity:
 *   - build: O(n log n) time and memory
 *   - lca: O(log n)
 *   - run: O(n sqrt(q) + q) calls to add and remove, plus O(q log q) sorting and O(q log n) for the lcas
 *
 * TreePathMo tree(n): nodes 0..n - 1 and no edges yet
 * tree.add_edge(u, v): adds the tree edge between u and v
 * tree.build(root = 0): roots the tree at root once all n - 1 edges are added
 * tree.lca(u, v): the lowest common ancestor of u and v, after build
 * tree.run(paths, add, remove, get, mode = NODES): res[i] = get(i) while the window holds the path paths[i] = {u, v}
 *
 * Notes:
 *   - add(node) and remove(node) receive tree nodes, get(i) returns the answer for the current window
 *   - mode = NODES: the window holds every node on the path u..v
 *   - mode = EDGES: the window holds the path without the lca, so store the value of edge (parent[v], v) on v
 *   - The window starts and ends empty, so the same callbacks and state can serve several runs
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY mo from code_library/data_structures/mo.cpp
template <typename Add, typename Remove, typename Get>
auto mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Get get){
    int q = queries.size();
    vector<decay_t<decltype(get(0))>> res(q);
    if (q == 0) return res;
    int block = max(1, (int)(n / sqrt((double)q)));

    /// Odd-even block order: Hilbert curve order measured about 5% more add / remove calls on random n = q = 5e5
    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? queries[a].second > queries[b].second : queries[a].second < queries[b].second;
    });

    /// Growing before shrinking keeps every removed position inside the window, empty targets included
    int cur_l = 0, cur_r = -1;
    for (int qi : order){
        auto [l, r] = queries[qi];
        assert(0 <= l && l <= r + 1 && r < n);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        res[qi] = get(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);
    return res;
}
/// END COPY mo

struct TreePathMo{
    enum Mode{NODES, EDGES};

    int n, lg;
    vector<vector<int>> adj, up;
    vector<int> depth, first, last, tour;

    TreePathMo(int n) : n(n), adj(n), depth(n), first(n), last(n){
        lg = 1;
        while ((1 << lg) < n) lg++;
    }

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    /// Euler tour with every node at its entry and exit, iterative so deep trees are safe
    void build(int root = 0){
        up.assign(lg + 1, vector<int>(n, root));
        tour.clear();
        vector<int> it(n, 0), stack = {root};
        depth[root] = 0, first[root] = 0, tour.push_back(root);

        while (!stack.empty()){
            int u = stack.back();
            if (it[u] < (int)adj[u].size()){
                int v = adj[u][it[u]++];
                if (v == up[0][u]) continue;
                up[0][v] = u, depth[v] = depth[u] + 1;
                first[v] = tour.size(), tour.push_back(v);
                stack.push_back(v);
                continue;
            }
            last[u] = tour.size(), tour.push_back(u);
            stack.pop_back();
        }

        assert((int)tour.size() == 2 * n);
        for (int k = 1; k <= lg; k++){
            for (int v = 0; v < n; v++) up[k][v] = up[k - 1][up[k - 1][v]];
        }
    }

    int lca(int u, int v) const{
        if (depth[u] < depth[v]) swap(u, v);
        for (int k = lg; k >= 0; k--){
            if (depth[u] - (1 << k) >= depth[v]) u = up[k][u];
        }

        if (u == v) return u;
        for (int k = lg; k >= 0; k--){
            if (up[k][u] != up[k][v]) u = up[k][u], v = up[k][v];
        }

        return up[0][u];
    }

    template <typename Add, typename Remove, typename Get>
    auto run(const vector<pair<int, int>>& paths, Add add, Remove remove, Get get, Mode mode = NODES){
        /// Tour positions toggle their node: [first[u], first[v]] if u is an ancestor, else [last[u], first[v]] + lca
        vector<pair<int, int>> ranges;
        vector<int> lca_to_add;  /// toggled in around get(i) for node paths whose range misses the lca, else -1
        for (auto [u, v] : paths){
            if (first[u] > first[v]) swap(u, v);
            int w = lca(u, v);
            if (w == u){
                /// Edge mode drops the lca u, u == v leaves it empty
                ranges.push_back({first[u] + (mode == EDGES), first[v]});
                lca_to_add.push_back(-1);
            }
            else{
                ranges.push_back({last[u], first[v]});
                lca_to_add.push_back(mode == EDGES ? -1 : w);
            }
        }

        vector<char> inside(n, 0);
        auto toggle = [&](int node){
            if (inside[node]) remove(node);
            else add(node);
            inside[node] ^= 1;
        };

        return mo(2 * n, ranges, [&](int i){ toggle(tour[i]); }, [&](int i){ toggle(tour[i]); }, [&](int i){
            if (lca_to_add[i] != -1) toggle(lca_to_add[i]);
            auto res = get(i);
            if (lca_to_add[i] != -1) toggle(lca_to_add[i]);
            return res;
        });
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
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
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
        for (; isdigit(c); c = read_char()){
            v = v * 10 - (c - '0');
        }
        x = neg ? v : -v;

        return true;
    }

    bool read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        if (!read_token()) return false;
        x = strtof(tmpbuf, nullptr);
        return true;
    }

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
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

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
        return true;
    }

    bool read_one(long double &x){
        if (!read_token()) return false;
        x = strtold(tmpbuf, nullptr);
        return true;
    }

    bool read_one(string &s){
        int c = ' ';

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !isspace(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (isspace(c));
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

#define MAX 100010
#define BLOCK 256

int w[MAX], cnt[MAX], present[MAX / BLOCK + 1];

void add_node(int v){
    if (w[v] < MAX && cnt[w[v]]++ == 0) present[w[v] / BLOCK]++;
}

void remove_node(int v){
    if (w[v] < MAX && --cnt[w[v]] == 0) present[w[v] / BLOCK]--;
}

int mex(int){
    int b = 0, y;
    while (present[b] == BLOCK) b++;
    y = b * BLOCK;
    while (cnt[y]) y++;

    return y;
}

int main(){
    int n, q, i, u, v;

    read(n, q);
    TreePathMo tree(n);
    vector<array<int, 3>> edges(n - 1);
    for (i = 0; i < n - 1; i++){
        read(edges[i][0], edges[i][1], edges[i][2]);
        tree.add_edge(--edges[i][0], --edges[i][1]);
    }
    tree.build();

    for (i = 0; i < n - 1; i++){
        u = edges[i][0], v = edges[i][1];
        w[tree.depth[u] > tree.depth[v] ? u : v] = edges[i][2];
    }

    vector<pair<int, int>> paths(q);
    for (i = 0; i < q; i++){
        read(u, v);
        paths[i] = {u - 1, v - 1};
    }

    auto res = tree.run(paths, add_node, remove_node, mex, TreePathMo::EDGES);
    for (i = 0; i < q; i++) write(res[i]);
    return 0;
}
