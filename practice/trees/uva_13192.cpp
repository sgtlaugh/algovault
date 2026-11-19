/***
 *
 * Lowest Common Ancestor
 * LCA, distance, jumps and path folds on a weighted tree, 0-based nodes
 *
 * Complexity:
 *   - LCA: O(n log n) to build, O(log n) per lca, dist, kth_ancestor and jump, O(n log n) memory
 *   - FastLCA: O(n log n) to build, O(1) per lca and dist, O(n log n) memory
 *   - PathFold: O(n log n) to build, O(log n) per query, n log n long longs of memory
 *   - path_intersection: 6 lca calls, O(1) with FastLCA, O(log n) with LCA
 *
 * LCA tree(n) or FastLCA tree(n): an empty tree on nodes 0..n-1, both share the calls below
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * tree.build(root): roots the tree at root, 0 by default, asserts n - 1 edges and every node reached exactly once
 * tree.depth(v): the number of edges from the root to v
 * tree.dist(u, v): the sum of the edge weights on the path u - v
 * tree.lca(u, v): the lowest common ancestor of u and v
 * tree.jump(u, v, k), LCA only: the k-th vertex on the path u -> v, u for k = 0, v for the last, -1 beyond v
 * tree.kth_ancestor(v, k), LCA only: the k-th ancestor of v, v itself for k = 0, -1 if k > depth(v)
 * tree.tin[v], FastLCA only: v's position in preorder, so v's subtree is a contiguous range starting at tin[v]
 * tree.order[i], FastLCA only: the vertex at preorder position i, so a reverse scan visits children before parents
 *
 * PathFold fold(tree, op): on a built LCA, op associative and commutative on long long, like max, min or gcd
 * fold.query(u, v): op over the edge weights on the path u - v, u != v
 *
 * path_intersection(tree, a, b, c, d): the common part of paths a - b and c - d as an unordered pair of endpoints
 *
 * Notes:
 *   - path_intersection takes either structure, returns {-1, -1} for no shared vertex, {x, x} for exactly x
 *   - Use FastLCA by default, its lca is 6x faster on random trees, 20x on paths; LCA for kth_ancestor, jump, PathFold
 *   - PathFold keeps a reference to tree, so build a new fold after tree.build with another root
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LCA{
    int n, lg, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<vector<int>> up;
    vector<int> level;
    vector<long long> weight_sum;

    LCA(int n) : n(n), adj(n), level(n), weight_sum(n){
        lg = 1;
        while ((1 << lg) < n) lg++;
    }

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    /// depth < 2^lg, so levels 0..lg-1 cover every jump
    void build(int root = 0){
        assert(added_edges == n - 1);
        up.assign(lg, vector<int>(n, -1));
        vector<int> order = {root};
        up[0][root] = root, level[root] = 0, weight_sum[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (auto e : adj[u]){ int v = e.first; long long w = e.second;
                if (v == up[0][u]) continue;
                assert(up[0][v] == -1);
                up[0][v] = u, level[v] = level[u] + 1, weight_sum[v] = weight_sum[u] + w;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int k = 1; k < lg; k++){
            for (int v = 0; v < n; v++) up[k][v] = up[k - 1][up[k - 1][v]];
        }
    }

    int depth(int v) const{
        return level[v];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }

    int lca(int u, int v) const{
        if (level[u] < level[v]) swap(u, v);
        u = kth_ancestor(u, level[u] - level[v]);
        if (u == v) return u;

        /// u != v at equal depth, so level[u] >= 1, and jumps longer than level[u] reach the root from both
        for (int k = __lg(level[u]); k >= 0; k--){
            if (up[k][u] != up[k][v]) u = up[k][u], v = up[k][v];
        }

        return up[0][u];
    }

    int jump(int u, int v, int k) const{
        int l = lca(u, v), up_len = level[u] - level[l], len = up_len + level[v] - level[l];
        if (k > len) return -1;
        return k <= up_len ? kth_ancestor(u, k) : kth_ancestor(v, len - k);
    }

    int kth_ancestor(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        for (int i = 0; k; i++, k >>= 1){
            if (k & 1) v = up[i][v];
        }

        return v;
    }
};

template <typename Op>
struct PathFold{
    const LCA& tree;
    Op op;
    vector<vector<long long>> table;

    /// table[k][v] folds the 2^k edges above v; jumps past the root fold in the root's 0 and are never queried
    PathFold(const LCA& tree, Op op) : tree(tree), op(op), table(tree.lg, vector<long long>(tree.n)){
        for (int v = 0; v < tree.n; v++) table[0][v] = tree.weight_sum[v] - tree.weight_sum[tree.up[0][v]];
        for (int k = 1; k < tree.lg; k++){
            for (int v = 0; v < tree.n; v++) table[k][v] = op(table[k - 1][v], table[k - 1][tree.up[k - 1][v]]);
        }
    }

    long long query(int u, int v) const{
        assert(u != v);
        int l = tree.lca(u, v), up_u = tree.level[u] - tree.level[l], up_v = tree.level[v] - tree.level[l];
        if (up_u == 0) return climb(v, up_v);
        if (up_v == 0) return climb(u, up_u);
        return op(climb(u, up_u), climb(v, up_v));
    }

    /// op over the k >= 1 edges above v
    long long climb(int v, int k) const{
        int i = __builtin_ctz(k);
        long long res = table[i][v];

        v = tree.up[i][v];
        for (k &= k - 1; k; k &= k - 1){
            i = __builtin_ctz(k);
            res = op(res, table[i][v]);
            v = tree.up[i][v];
        }

        return res;
    }
};

/// BEGIN SHARED fast_lca
struct FastLCA{
    int n, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<int> tin, order, level;
    vector<long long> weight_sum;
    vector<vector<int>> table;

    FastLCA(int n) : n(n), adj(n), tin(n), order(n), level(n), weight_sum(n) {}

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    /// table[k][i] is the min tin over the parents of the vertices at preorder positions i..i+2^k-1
    void build(int root = 0){
        assert(added_edges == n - 1);
        vector<int> parent(n, -1);
        vector<pair<int, int>> stack = {{root, 0}};
        parent[root] = root, level[root] = 0, weight_sum[root] = 0;
        table.assign(1, vector<int>(n));
        int visited = 0;

        while (!stack.empty()){
            int u = stack.back().first, parent_tin = stack.back().second;
            stack.pop_back();
            tin[u] = visited, order[visited] = u, table[0][visited++] = parent_tin;
            for (auto e : adj[u]){ int v = e.first; long long w = e.second;
                if (v == parent[u]) continue;
                assert(parent[v] == -1);
                parent[v] = u, level[v] = level[u] + 1, weight_sum[v] = weight_sum[u] + w;
                stack.push_back({v, tin[u]});
            }
        }
        assert(visited == n);

        for (int k = 1; (1 << k) <= n; k++){
            table.push_back(vector<int>(n - (1 << k) + 1));
            for (int i = 0; i + (1 << k) <= n; i++){
                table[k][i] = min(table[k - 1][i], table[k - 1][i + (1 << (k - 1))]);
            }
        }
    }

    int depth(int v) const{
        return level[v];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }

    /// For tin[u] < tin[v], the shallowest vertex at preorder positions tin[u] + 1..tin[v] is a child of the lca
    int lca(int u, int v) const{
        if (u == v) return u;
        int l = tin[u], r = tin[v];
        if (l > r) swap(l, r);
        int k = __lg(r - l);
        return order[min(table[k][l + 1], table[k][r - (1 << k) + 1])];
    }
};
/// END SHARED fast_lca

/// The two deepest cross lcas are the endpoints, the paths meet iff the deepest is not above the deeper path top
template <typename Tree>
pair<int, int> path_intersection(const Tree& tree, int a, int b, int c, int d){
    int top = tree.lca(a, b), other_top = tree.lca(c, d);
    if (tree.depth(other_top) > tree.depth(top)) top = other_top;

    array<int, 4> cand = {tree.lca(a, c), tree.lca(a, d), tree.lca(b, c), tree.lca(b, d)};
    sort(cand.begin(), cand.end(), [&](int x, int y){ return tree.depth(x) > tree.depth(y); });
    if (tree.depth(cand[0]) < tree.depth(top)) return {-1, -1};
    return {cand[0], cand[1]};
}

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

#define MAX 50010
#define V 101

int n, A[MAX], B[MAX], q[MAX], parent[MAX], pedge[MAX];
double L[MAX][V];
vector<pair<int, int>> adj[MAX];

struct Max{
    long long operator()(long long x, long long y) const{
        return max(x, y);
    }
};

void prefix(){
    int i, v, head = 0, tail = 0;

    for (i = 0; i < n; i++) parent[i] = -1;
    for (v = 0; v < V; v++) L[0][v] = 0;
    q[tail++] = 0, parent[0] = 0;
    while (head < tail){
        int u = q[head++];
        for (i = 0; i < (int)adj[u].size(); i++){
            int x = adj[u][i].first, e = adj[u][i].second;
            if (parent[x] != -1) continue;
            parent[x] = u, q[tail++] = x;
            for (v = 0; v < V; v++){
                double cur = 0;
                if (v > A[e] && v <= B[e]) cur = log((double)(v - A[e]) / (B[e] - A[e] + 1));
                L[x][v] = L[u][v] + cur;
            }
        }
    }
}

int main(){
    int i, u, v, m, x, y;

    while (read(n)){
        LCA tree(n);
        for (i = 0; i < n; i++) adj[i].clear();
        for (i = 0; i < n - 1; i++){
            read(u, v, A[i], B[i]);
            u--, v--;
            tree.add_edge(u, v, A[i]);
            adj[u].push_back(make_pair(v, i));
            adj[v].push_back(make_pair(u, i));
        }
        tree.build(0);
        PathFold<Max> fold(tree, Max());
        prefix();

        read(m);
        while (m--){
            read(x, y);
            x--, y--;
            int l = tree.lca(x, y), top = fold.query(x, y);

            double res = top;
            for (v = top + 1; v < V; v++) res += 1.0 - exp(L[x][v] + L[y][v] - 2 * L[l][v]);
            write(fixed(res, 9));
        }
    }

    return 0;
}
