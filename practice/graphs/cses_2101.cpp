/***
 *
 * Kruskal Reconstruction Tree
 * The vertices reachable from v using only edges of weight <= w, as one contiguous range of a fixed vertex order
 *
 * Complexity: O(m log m + n log n) to build, O(log n) per query, O(n log n) memory
 *
 * KruskalReconstructionTree<T> krt(n): vertices 0 to n - 1, no edges
 * krt.add_edge(u, v, w): undirected edge of weight w, every edge goes in before build
 * krt.build(): builds the tree
 * krt.node(v, w): the tree node of v's component under edges <= w, v itself when it is alone
 * krt.positions(v, w): the range [l, r] of positions in krt.order holding exactly v's component
 * krt.size(v, w): the number of vertices in v's component
 * krt.weight[x]: for a merge node x >= n, the weight of the edge that created it, the component's bottleneck
 * krt.order[i]: the vertex at position i
 * krt.pos[v]: the position of vertex v
 *
 * Notes:
 *   - Vertices are leaves 0 to n - 1, merge nodes get ids n, n + 1, ..., up to 2n - 2
 *   - Disconnected graphs, self loops, multi-edges and negative weights are fine
 *   - A segment tree over krt.order aggregates everything reachable from v with edges <= w
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T = long long>
struct KruskalReconstructionTree{
    int n, lg = 1;
    vector<tuple<T, int, int>> edges;
    vector<T> weight;
    vector<int> order, pos, first, leaf_count;
    vector<vector<int>> up;

    KruskalReconstructionTree(int n) : n(n) {}

    void add_edge(int u, int v, T w){
        edges.push_back({w, u, v});
    }

    void build(){
        sort(edges.begin(), edges.end());
        vector<int> dsu(n), parent(n);
        vector<array<int, 2>> child(n, {-1, -1});
        iota(dsu.begin(), dsu.end(), 0);
        iota(parent.begin(), parent.end(), 0);
        weight.assign(n, T());

        auto find = [&](int x){
            int root = x;
            while (dsu[root] != root) root = dsu[root];
            while (dsu[x] != root){
                int next = dsu[x];
                dsu[x] = root, x = next;
            }

            return root;
        };

        for (auto [w, u, v] : edges){
            u = find(u), v = find(v);
            if (u == v) continue;
            int id = weight.size();
            dsu[u] = dsu[v] = parent[u] = parent[v] = id;
            dsu.push_back(id), parent.push_back(id), child.push_back({u, v}), weight.push_back(w);
        }
        int m = weight.size();

        /// children always have smaller ids than their parent, so id order is bottom-up and reverse id order top-down
        leaf_count.assign(m, 1);
        for (int x = n; x < m; x++) leaf_count[x] = leaf_count[child[x][0]] + leaf_count[child[x][1]];

        first.assign(m, -1);
        int offset = 0;
        for (int x = m - 1; x >= 0; x--){
            if (parent[x] == x) first[x] = offset, offset += leaf_count[x];
            if (x >= n) first[child[x][0]] = first[x], first[child[x][1]] = first[x] + leaf_count[child[x][0]];
        }

        order.assign(n, 0), pos.assign(n, 0);
        for (int v = 0; v < n; v++) pos[v] = first[v], order[first[v]] = v;

        while ((1 << lg) < n) lg++;
        up.assign(lg, parent);
        for (int k = 1; k < lg; k++){
            for (int x = 0; x < m; x++) up[k][x] = up[k - 1][up[k - 1][x]];
        }
    }

    /// a root is its own parent, so jumping past the root lands on the root, which is correct when its weight is <= w
    int node(int v, T w) const{
        for (int k = lg - 1; k >= 0; k--){
            if (weight[up[k][v]] <= w) v = up[k][v];
        }

        return v;
    }

    pair<int, int> positions(int v, T w) const{
        int x = node(v, w);
        return {first[x], first[x] + leaf_count[x] - 1};
    }

    int size(int v, T w) const{
        return leaf_count[node(v, w)];
    }
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
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
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
            if (buf_len == 0) return EOF;
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

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

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

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
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

KruskalReconstructionTree<int> krt(0);
vector<int> depth;

int lca(int a, int b){
    if (depth[a] < depth[b]) swap(a, b);
    for (int k = krt.lg - 1; k >= 0; k--){
        if ((depth[a] - (1 << k)) >= depth[b]) a = krt.up[k][a];
    }
    if (a == b) return a;
    for (int k = krt.lg - 1; k >= 0; k--){
        if (krt.up[k][a] != krt.up[k][b]) a = krt.up[k][a], b = krt.up[k][b];
    }

    return (krt.up[0][a] == krt.up[0][b]) ? krt.up[0][a] : -1;
}

int main(){
    int n, m, q, i, a, b;

    read(n, m, q);
    krt = KruskalReconstructionTree<int>(n);
    for (i = 1; i <= m; i++){
        read(a, b);
        krt.add_edge(a - 1, b - 1, i);
    }
    krt.build();

    int k = krt.weight.size();
    depth.assign(k, 0);
    for (i = k - 1; i >= 0; i--){
        if (krt.up[0][i] != i) depth[i] = depth[krt.up[0][i]] + 1;
    }

    while (q--){
        read(a, b);
        a--, b--;
        if (a == b){
            write(0);
            continue;
        }

        int c = lca(a, b);
        write(c == -1 ? -1 : krt.weight[c]);
    }

    return 0;
}
