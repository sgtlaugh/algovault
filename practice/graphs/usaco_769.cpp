/***
 *
 * Articulation Points, Biconnected Components and the Block-Cut Tree
 * Splits an undirected graph into its blocks (vertex-biconnected components), parallel edges and self loops allowed
 *
 * Complexity: O(n + m)
 *
 * Biconnected g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v): adds an edge and returns its id, ids count up from 0 in add_edge order
 * g.run(): fills blocks, cuts and is_cut, rerun it after more add_edge calls
 * g.block_cut_tree(): adjacency lists over n + blocks nodes, node v < n is vertex v, node n + b is block b, after run()
 * g.blocks: edge ids of each block, two edges share a block iff they lie on a common simple cycle
 * g.cuts: sorted articulation points, the vertices whose removal increases the number of connected components
 * g.is_cut[v]: whether v is an articulation point
 *
 * Notes:
 *   - Every block links to each vertex it touches, so cut vertices have degree >= 2 and the rest are leaves or isolated
 *   - The block-cut tree is a forest with one tree per connected component
 *   - The vertices separating u from w are exactly the cut vertices strictly inside the tree path from u to w
 *   - A bridge is a block of its own, a doubled bridge is a 2-edge block, a self loop belongs to no block
 *   - A vertex lies in several blocks iff it is an articulation point
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Biconnected{
    int n;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge id)
    vector<pair<int, int>> edges;
    vector<vector<int>> blocks;
    vector<int> cuts;
    vector<bool> is_cut;

    Biconnected(int n) : n(n), adj(n) {}

    int add_edge(int u, int v){
        int id = edges.size();
        edges.push_back({u, v});
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        return id;
    }

    void run(){
        vector<int> discover(n, 0), low(n), parent_edge(n, -1), next_edge(n, 0), path, edge_stack;
        is_cut.assign(n, false);
        blocks.clear(), cuts.clear();

        int disc_t = 0;
        for (int root = 0; root < n; root++){
            if (discover[root]) continue;
            int children = 0;
            discover[root] = low[root] = ++disc_t;
            path.push_back(root);

            while (!path.empty()){
                int u = path.back();
                if (next_edge[u] < (int)adj[u].size()){
                    auto [v, id] = adj[u][next_edge[u]++];
                    /// Skip only the edge we arrived by, a parallel edge back to the parent is a real back edge
                    if (v == u || id == parent_edge[u]) continue;
                    if (!discover[v]){
                        discover[v] = low[v] = ++disc_t;
                        parent_edge[v] = id;
                        edge_stack.push_back(id);
                        path.push_back(v);
                    }
                    else if (discover[v] < discover[u]){
                        low[u] = min(low[u], discover[v]);
                        edge_stack.push_back(id);
                    }
                    continue;
                }

                path.pop_back();
                if (path.empty()) break;
                int p = path.back();
                low[p] = min(low[p], low[u]);
                if (low[u] < discover[p]) continue;

                if (p != root) is_cut[p] = true;
                else if (++children == 2) is_cut[p] = true;

                blocks.emplace_back();
                int id;
                do{
                    id = edge_stack.back();
                    edge_stack.pop_back();
                    blocks.back().push_back(id);
                } while (id != parent_edge[u]);
            }
        }

        for (int u = 0; u < n; u++){
            if (is_cut[u]) cuts.push_back(u);
        }
    }

    vector<vector<int>> block_cut_tree(){
        vector<vector<int>> tree(n + blocks.size());
        vector<int> last(n, -1);
        for (int b = 0; b < (int)blocks.size(); b++){
            for (int id: blocks[b]){
                for (int u: {edges[id].first, edges[id].second}){
                    if (last[u] == b) continue;
                    last[u] = b;
                    tree[u].push_back(n + b);
                    tree[n + b].push_back(u);
                }
            }
        }

        return tree;
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

#define MAX 1510

const int dx[] = {-1, 0, 1, 0}, dy[] = {0, 1, 0, -1};

int n, m;
string grid[MAX];
vector<int> eid[4], block, cow;
vector<char> seen, good;

bool valid(int x, int y){
    return x >= 0 && x < n && y >= 0 && y < m && grid[x][y] != '#';
}

void walk(int sx, int sy, int bx, int by){
    vector<int> q;

    cow.assign(n * m, 0);
    cow[sx * m + sy] = 1, q.push_back(sx * m + sy);
    for (int i = 0; i < (int)q.size(); i++){
        int x = q[i] / m, y = q[i] % m;
        for (int d = 0; d < 4; d++){
            int u = x + dx[d], v = y + dy[d];
            if (!valid(u, v) || (u == bx && v == by) || cow[u * m + v]) continue;
            cow[u * m + v] = 1, q.push_back(u * m + v);
        }
    }
}

int main(){
    int q, i, j, d, e, x, y, ax = 0, ay = 0, bx = 0, by = 0;

    freopen("pushabox.in", "r", stdin);
    freopen("pushabox.out", "w", stdout);

    read(n, m, q);
    for (i = 0; i < n; i++){
        read(grid[i]);
        for (j = 0; j < m; j++){
            if (grid[i][j] == 'A') ax = i, ay = j;
            if (grid[i][j] == 'B') bx = i, by = j;
        }
    }

    Biconnected g(n * m);
    for (d = 0; d < 4; d++) eid[d].assign(n * m, -1);
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            if (!valid(i, j)) continue;
            for (d = 1; d <= 2; d++){
                x = i + dx[d], y = j + dy[d];
                if (!valid(x, y)) continue;
                e = g.add_edge(i * m + j, x * m + y);
                eid[d][i * m + j] = e, eid[(d + 2) & 3][x * m + y] = e;
            }
        }
    }

    g.run();
    block.assign(g.edges.size(), -1);
    for (i = 0; i < (int)g.blocks.size(); i++){
        for (int id: g.blocks[i]) block[id] = i;
    }

    walk(ax, ay, bx, by);
    seen.assign(4 * n * m, 0), good.assign(n * m, 0);
    good[bx * m + by] = 1;

    vector<int> st;
    for (d = 0; d < 4; d++){
        x = bx + dx[d], y = by + dy[d];
        if (valid(x, y) && cow[x * m + y]) seen[(bx * m + by) * 4 + d] = 1, st.push_back((bx * m + by) * 4 + d);
    }

    for (i = 0; i < (int)st.size(); i++){
        int c = st[i] >> 2, k = st[i] & 3;
        x = c / m, y = c % m;
        good[c] = 1;

        for (d = 0; d < 4; d++){
            if (d == k || eid[d][c] == -1 || block[eid[d][c]] != block[eid[k][c]]) continue;
            if (!seen[c * 4 + d]) seen[c * 4 + d] = 1, st.push_back(c * 4 + d);
        }

        int u = x - dx[k], v = y - dy[k];
        if (valid(u, v) && !seen[(u * m + v) * 4 + k]) seen[(u * m + v) * 4 + k] = 1, st.push_back((u * m + v) * 4 + k);
    }

    while (q--){
        read(x, y);
        write(good[(x - 1) * m + (y - 1)] ? "YES" : "NO");
    }

    return 0;
}
