/***
 *
 * Static Top Tree
 * A tree DP kept up to date under point updates, answered from the fixed root or from any vertex as the root
 *
 * Complexity: O(n log n) time and O(n) operation calls to build, O(log n) operation calls per update and per value
 *
 * StaticTopTree<Ops> tree(n, edges, ops, root = 0, roots = ANY_ROOT): vertices 0..n-1, edge e joins edges[e]
 * tree.update(v): recomputes after the data of vertex v changed
 * tree.update_edge(e): recomputes after the data of edge e changed
 * tree.value(): the DP of the whole tree rooted at root, a Path
 * tree.value(r): the DP of the whole tree rooted at r, a Path built as ops.add_vertex(..., r, -1) or vertex(r, -1)
 *
 * Notes:
 *   - A Path is a chain of vertices x_1..x_k, everything hanging off them, and the edge e joining x_1 to the outside
 *   - A Path's last vertex x_k may join the outside by an edge the Path does not hold, and e = -1 means none
 *   - A Point is everything hanging off one vertex through one edge, that edge included
 *   - Ops defines types Path and Point, both default constructible, and five calls, read from the data they hold:
 *       - vertex(v, e): the Path of v alone with edge e
 *       - add_vertex(x, v, e): the Path of v with edge e and the Point x hanging off v
 *       - add_edge(p): p as a Point hanging off the vertex its edge joins, p's last vertex joining nothing beyond
 *       - rake(a, b): two Points hanging off the same vertex as one, associative and commutative
 *       - compress(top, bottom): one Path, the edge of bottom joining its first vertex to the last vertex of top
 *   - value(r) joins Paths built from both ends of each chain, so the five calls must describe one consistent DP
 *   - ANY_ROOT keeps every Path in both directions; FIXED_ROOT keeps one, 15-30% faster, and value(r) then asserts
 *   - Ops is copied, so it should reach mutable data through pointers or references
 *   - Recursion depth is O(log n), deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Ops>
struct StaticTopTree{
    using Path = typename Ops::Path;
    using Point = typename Ops::Point;
    enum Kind{VERTEX, ADD_VERTEX, COMPRESS, ADD_EDGE, RAKE};
    enum Roots{ANY_ROOT, FIXED_ROOT};

    int n, root, top_node;
    Roots roots;
    Ops ops;
    vector<pair<int, int>> edges;
    vector<int> parent, parent_edge, heavy, heavy_edge, size, top, vertex_node;
    vector<char> kind;
    vector<int> left, right, up, vtx;
    vector<Path> fwd, rev;  /// rev is the same Path seen from its other end, edge heavy_edge instead of parent_edge
    vector<Point> point;

    StaticTopTree(int n, const vector<pair<int, int>>& edges, Ops ops, int root = 0, Roots roots = ANY_ROOT)
        : n(n), root(root), roots(roots), ops(ops), edges(edges), parent(n, -1), parent_edge(n, -1), heavy(n, -1), heavy_edge(n, -1),
          size(n, 1), top(n), vertex_node(n){
        assert(n >= 1 && (int)edges.size() == n - 1 && 0 <= root && root < n);

        vector<int> start(n + 1, 0), adj(2 * (n - 1));
        for (auto [u, v] : edges){
            assert(0 <= u && u < n && 0 <= v && v < n);
            start[u + 1]++, start[v + 1]++;
        }
        for (int v = 0; v < n; v++) start[v + 1] += start[v];
        vector<int> pos(start.begin(), start.end() - 1);
        for (int e = 0; e < n - 1; e++) adj[pos[edges[e].first]++] = e, adj[pos[edges[e].second]++] = e;

        vector<int> order = {root};
        vector<char> seen(n, 0);
        seen[root] = 1;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int e = adj[j], v = edges[e].first ^ edges[e].second ^ u;
                if (seen[v]) continue;
                seen[v] = 1, parent[v] = u, parent_edge[v] = e;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i > 0; i--) size[parent[order[i]]] += size[order[i]];
        for (int i = 1; i < n; i++){
            int v = order[i], p = parent[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v, heavy_edge[p] = parent_edge[v];
        }

        vector<vector<int>> light(n);
        top[root] = root;
        for (int i = 1; i < n; i++){
            int v = order[i];
            top[v] = heavy[parent[v]] == v ? top[parent[v]] : v;
            if (heavy[parent[v]] != v) light[parent[v]].push_back(v);
        }
        top_node = build_chain(root, light);
        up[top_node] = -1;
    }

    void update(int v){
        for (int node = vertex_node[v]; node != -1; node = up[node]) recompute(node);
    }

    /// Only the child's Path holds e, plus the parent's reversed Path when e is its heavy edge
    void update_edge(int e){
        auto [u, v] = edges[e];
        int child = parent[u] == v ? u : v;
        update(child);
        if (roots == ANY_ROOT && heavy_edge[parent[child]] == e) update(parent[child]);
    }

    Path value() const{
        return fwd[top_node];
    }

    /// Climbs from r: each chain on the way gives the part below and above the vertex the climb entered it at
    Path value(int r){
        assert(roots == ANY_ROOT);

        struct Level{
            int v;
            bool has_light = false, has_below = false, has_above = false;
            Point light;
            Path below, above;  /// below is a forward Path, above a reversed one starting next to v
        };

        vector<Level> levels;
        Level cur;
        cur.v = r;
        int node = vertex_node[r];
        if (kind[node] == ADD_VERTEX) cur.has_light = true, cur.light = point[left[node]];

        for (int child = node, p = up[node]; p != -1; child = p, p = up[p]){
            if (kind[p] == COMPRESS && child == left[p]){
                cur.below = cur.has_below ? ops.compress(cur.below, fwd[right[p]]) : fwd[right[p]];
                cur.has_below = true;
            }
            else if (kind[p] == COMPRESS){
                cur.above = cur.has_above ? ops.compress(cur.above, rev[left[p]]) : rev[left[p]];
                cur.has_above = true;
            }
            else if (kind[p] == ADD_EDGE){
                levels.push_back(cur);
                cur = Level();
            }
            else if (kind[p] == RAKE){
                const Point& other = point[child == left[p] ? right[p] : left[p]];
                cur.light = cur.has_light ? ops.rake(cur.light, other) : other;
                cur.has_light = true;
            }
            else cur.v = vtx[p];
        }
        levels.push_back(cur);

        /// rest is everything beyond the top of the chain of levels[i], a Path whose edge joins that chain's top
        bool has_rest = levels.back().has_above;
        Path rest = levels.back().above;
        for (int i = (int)levels.size() - 2; i >= -1; i--){
            const Level& at = levels[i + 1];
            bool has_x = at.has_light;
            Point x = at.light;
            auto hang = [&](const Path& p){
                Point hanging = ops.add_edge(p);
                x = has_x ? ops.rake(x, hanging) : hanging;
                has_x = true;
            };
            if (at.has_below) hang(at.below);
            if (has_rest) hang(rest);

            if (i == -1) return has_x ? ops.add_vertex(x, at.v, -1) : ops.vertex(at.v, -1);

            int e = parent_edge[top[levels[i].v]];
            Path beyond = has_x ? ops.add_vertex(x, at.v, e) : ops.vertex(at.v, e);
            rest = levels[i].has_above ? ops.compress(levels[i].above, beyond) : beyond;
            has_rest = true;
        }

        return rest;
    }

    int build_chain(int head, const vector<vector<int>>& light){
        vector<int> nodes, weights;
        for (int v = head; v != -1; v = heavy[v]){
            vector<int> points, point_weights;
            for (int c : light[v]){
                points.push_back(add_node(ADD_EDGE, build_chain(c, light), -1, -1));
                point_weights.push_back(size[c]);
            }

            int node;
            if (points.empty()) node = add_node(VERTEX, -1, -1, v);
            else node = add_node(ADD_VERTEX, merge(points, point_weights, 0, points.size(), RAKE), -1, v);
            vertex_node[v] = node;
            nodes.push_back(node), weights.push_back(size[v] - (heavy[v] == -1 ? 0 : size[heavy[v]]));
        }

        return merge(nodes, weights, 0, nodes.size(), COMPRESS);
    }

    /// Splits at the weight midpoint, so a vertex lies under O(log n) nodes in all
    int merge(const vector<int>& nodes, const vector<int>& weights, int l, int r, Kind k){
        if (r - l == 1) return nodes[l];

        long long total = 0, prefix = 0;
        for (int i = l; i < r; i++) total += weights[i];
        int m = l + 1;
        for (prefix = weights[l]; m < r - 1 && 2 * (prefix + weights[m]) <= total; m++) prefix += weights[m];

        return add_node(k, merge(nodes, weights, l, m, k), merge(nodes, weights, m, r, k), -1);
    }

    int add_node(Kind k, int l, int r, int v){
        int node = kind.size();
        kind.push_back(k), left.push_back(l), right.push_back(r), up.push_back(-1), vtx.push_back(v);
        fwd.emplace_back(), rev.emplace_back(), point.emplace_back();
        if (l != -1) up[l] = node;
        if (r != -1) up[r] = node;
        recompute(node);
        return node;
    }

    void recompute(int node){
        int l = left[node], r = right[node], v = vtx[node];
        bool both = roots == ANY_ROOT;
        if (kind[node] == VERTEX){
            fwd[node] = ops.vertex(v, parent_edge[v]);
            if (both) rev[node] = ops.vertex(v, heavy_edge[v]);
        }
        else if (kind[node] == ADD_VERTEX){
            fwd[node] = ops.add_vertex(point[l], v, parent_edge[v]);
            if (both) rev[node] = ops.add_vertex(point[l], v, heavy_edge[v]);
        }
        else if (kind[node] == COMPRESS){
            fwd[node] = ops.compress(fwd[l], fwd[r]);
            if (both) rev[node] = ops.compress(rev[r], rev[l]);
        }
        else if (kind[node] == ADD_EDGE) point[node] = ops.add_edge(fwd[l]);
        else point[node] = ops.rake(point[l], point[r]);
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

#define MAX 200010

const long long MOD = 998244353;

int n, children[MAX];
long long val[MAX];

struct Hash{
    struct Path{ long long a = 1, b = 0; };
    struct Point{ long long x = 1; };

    Path vertex(int v, int) const{ return {children[v] ? 1 : 0, val[v]}; }
    Path add_vertex(const Point& p, int v, int) const{ return {p.x, val[v]}; }
    Point add_edge(const Path& p) const{ return {p.b}; }
    Point rake(const Point& p, const Point& q) const{ return {p.x * q.x % MOD}; }
    Path compress(const Path& top, const Path& bottom) const{
        return {top.a * bottom.a % MOD, (top.a * bottom.b + top.b) % MOD};
    }
};

int main(){
    int q, i, p, v;

    read(n, q);
    vector<pair<int, int>> edges;
    for (i = 1; i < n; i++){
        read(p);
        edges.push_back({p - 1, i});
        children[p - 1]++;
    }
    for (i = 0; i < n; i++) read(val[i]);

    StaticTopTree<Hash> tree(n, edges, Hash(), 0, StaticTopTree<Hash>::FIXED_ROOT);
    while (q--){
        read(v);
        read(val[v - 1]);
        tree.update(v - 1);
        write(tree.value().b);
    }

    return 0;
}
