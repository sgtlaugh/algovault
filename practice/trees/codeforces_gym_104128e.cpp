/***
 *
 * Virtual Tree (Auxiliary Tree)
 * Compresses a rooted tree to a vertex subset plus the pairwise LCAs needed to keep its shape, 0-based nodes
 *
 * Complexity: O(n log n) to build and memory, O(k log k) per compress of k vertices
 *
 * VirtualTree tree(n): an empty tree on nodes 0..n-1
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * tree.build(root): roots the tree at root, 0 by default, asserts n - 1 edges and every node reached exactly once
 * tree.compress(nodes): the virtual tree of nodes as a list of VirtualTree::Node {v, parent, len}
 * tree.depth(v): the number of edges from the root to v
 * tree.dist(u, v): the sum of the edge weights on the path u - v
 * tree.lca(u, v): the lowest common ancestor of u and v
 *
 * Notes:
 *   - parent is the index of v's parent in the list, -1 for the root at index 0
 *   - len is the original length of the compressed edge, dist(v, list[parent].v), 0 for the root
 *   - The list is in DFS order, so parents come before children and every subtree is a contiguous range
 *   - The list has at most 2k - 1 vertices, duplicates in nodes are ignored, an empty nodes gives an empty list
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY fast_lca from code_library/trees/lca.cpp
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
            auto [u, parent_tin] = stack.back();
            stack.pop_back();
            tin[u] = visited, order[visited] = u, table[0][visited++] = parent_tin;
            for (auto [v, w] : adj[u]){
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
/// END COPY fast_lca

struct VirtualTree : FastLCA{
    using FastLCA::FastLCA;

    struct Node{ int v, parent; long long len; };

    /// The parent of each vertex is the lca of it and its predecessor in DFS order, so no stack is needed
    vector<Node> compress(vector<int> nodes) const{
        auto by_tin = [&](int u, int v){ return tin[u] < tin[v]; };

        sort(nodes.begin(), nodes.end(), by_tin);
        int k = nodes.size();
        for (int i = 0; i + 1 < k; i++) nodes.push_back(lca(nodes[i], nodes[i + 1]));
        sort(nodes.begin(), nodes.end(), by_tin);
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.empty()) return {};

        vector<Node> res = {{nodes[0], -1, 0}};
        for (int i = 1; i < (int)nodes.size(); i++){
            int p = lca(nodes[i - 1], nodes[i]);
            int p_idx = lower_bound(nodes.begin(), nodes.end(), p, by_tin) - nodes.begin();
            res.push_back({nodes[i], p_idx, weight_sum[nodes[i]] - weight_sum[p]});
        }

        return res;
    }
};

/***
 *
 * Sparse Table
 * O(1) range queries on a static array for operations that ignore repeats: min, max, gcd, and, or
 *
 * Complexity:
 *   - SparseTable: O(n log n) build time and memory, O(1) per query
 *   - LinearSparseTable: O(n) build time and memory, O(1) per query
 *
 * SparseTable<T> table(a): range minimum of a
 * SparseTable table(a, op): range op of a, for an op that ignores repeats like max or gcd
 * LinearSparseTable<T> table(a): range minimum of a
 * LinearSparseTable<T, greater<T>> table(a): range maximum of a
 * table.query(l, r): the answer over a[l..r], 0 <= l <= r < n
 *
 * Notes:
 *   - Sum, xor and product count repeats, use disjoint_sparse_table.cpp for them
 *   - At n = 1e6 LinearSparseTable builds 3x faster in 12 MB instead of 72 MB, and its queries are about 2x slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct MinOp{
    T operator()(const T& x, const T& y) const{
        return min(x, y);
    }
};

template <typename T, typename Op = MinOp<T>>
struct SparseTable{
    vector<vector<T>> dp;  /// dp[h][i] = op over a[i .. i + 2^h - 1]
    Op combine;

    SparseTable(const vector<T>& a, Op combine = Op()) : dp(1, a), combine(combine){
        int n = a.size();
        for (int h = 1; (1 << h) <= n; h++){
            dp.push_back(vector<T>(n - (1 << h) + 1));
            for (int i = 0; i + (1 << h) <= n; i++) dp[h][i] = combine(dp[h - 1][i], dp[h - 1][i + (1 << (h - 1))]);
        }
    }

    /// Two blocks of length 2^h cover [l, r], overlapping in the middle
    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)dp[0].size());
        int h = __lg(r - l + 1);
        return combine(dp[h][l], dp[h][r - (1 << h) + 1]);
    }
};

template <typename T, typename Compare = less<T>>
struct LinearSparseTable{
    static constexpr int B = 64;
    Compare compare;
    vector<T> val;
    vector<unsigned long long> mask;  /// bit k of mask[i]: val[i - k] beats every value after it up to i
    vector<vector<int>> table;        /// sparse table of block winners, by index

    LinearSparseTable(const vector<T>& a) : val(a), mask(a.size()){
        int n = val.size();
        unsigned long long stack = 0;  /// the monotonic stack of the last 64 positions, as a bitmask

        for (int i = 0; i < n; i++){
            stack <<= 1;
            /// the lowest set bit is the newest entry, pop while val[i] is at least as good
            while (stack && !compare(val[i - __builtin_ctzll(stack)], val[i])) stack &= stack - 1;
            stack |= 1;
            mask[i] = stack;
        }

        int blocks = n / B;  /// a query reads only blocks wholly before r's, so a partial last block is never read
        table.assign(1, vector<int>(blocks));
        for (int b = 0; b < blocks; b++) table[0][b] = best_ending_at(b * B + B - 1);
        for (int k = 1; (1 << k) <= blocks; k++){
            table.push_back(vector<int>(blocks - (1 << k) + 1));
            for (int b = 0; b + (1 << k) <= blocks; b++) table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
        }
    }

    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)val.size());
        if (r - l + 1 <= B) return val[best_ending_at(r, r - l + 1)];

        int res = better(best_ending_at(l + B - 1), best_ending_at(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = __lg(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }

        return val[res];
    }

    /// Index of the best value in [r - len + 1, r], len <= 64: the oldest stack entry inside the window
    int best_ending_at(int r, int len = B) const{
        unsigned long long m = len == B ? mask[r] : mask[r] & ((1ULL << len) - 1);
        return r - __lg(m);
    }

    int better(int i, int j) const{
        return compare(val[j], val[i]) ? j : i;
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
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - Floats print at printf's speed, a floating token read keeps its first 127 characters
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

    bool read_one(double &x){
        if (!read_token()) return false;
        x = strtod(tmpbuf, nullptr);
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

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
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

int n;
vector<int> groups[MAX];

long long solve(const VirtualTree& tree, const SparseTable<long long>& table, int d){
    int i;
    vector<VirtualTree::Node> vt = tree.compress(groups[d]);
    vector<long long> sum(vt.size(), 0), f(vt.size());

    for (i = (int)vt.size() - 1; i >= 0; i--){
        int v = vt[i].v, top = 0;
        if (vt[i].parent != -1) top = tree.depth(vt[vt[i].parent].v) + 1;

        f[i] = table.query(d - tree.depth(v), d - top);
        if (tree.depth(v) != d) f[i] = min(f[i], sum[i]);
        if (vt[i].parent != -1) sum[vt[i].parent] += f[i];
    }

    return f[0];
}

int main(){
    int t, i, u, v;

    read(t);
    while (t--){
        read(n);
        vector<long long> cost(n);
        read(cost);

        VirtualTree tree(n);
        for (i = 1; i < n; i++){
            read(u, v);
            tree.add_edge(u - 1, v - 1);
        }
        tree.build(0);

        SparseTable<long long> table(cost);
        for (i = 0; i < n; i++) groups[tree.depth(i)].push_back(i);

        long long res = 0;
        for (i = 0; i < n && !groups[i].empty(); i++){
            res += solve(tree, table, i);
            groups[i].clear();
        }
        write(res);
    }

    return 0;
}
