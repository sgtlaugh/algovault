/***
 *
 * Gomory-Hu Tree (Gusfield)
 * All-pairs minimum cut of an undirected graph from n - 1 max flows
 *
 * Complexity:
 *   - build: n - 1 Dinic max flows, each O(n^2 m), O(n + m) memory
 *   - min_cut: O(n)
 *
 * GomoryHu g(n): nodes 0 to n - 1, no edges
 * g.add_edge(u, v, cap): undirected edge of capacity cap
 * g.build(): builds the tree, call it after all edges are added
 * g.min_cut(u, v): the minimum cut between u != v
 * g.tree_parent[v]: for v >= 1, v's parent in the tree rooted at 0
 * g.tree_weight[v]: for v >= 1, the min cut between v and tree_parent[v]
 *
 * Notes:
 *   - Removing tree edge (v, tree_parent[v]) leaves the subtree of v as one side of a minimum cut between its endpoints
 *   - Parallel edges and self loops are allowed, a disconnected pair has min cut 0
 *   - Total capacity must fit in a long long
 *   - Dinic's dfs recurses up to n deep
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GomoryHu{
    struct Edge{
        int v;
        long long cap, flow;
    };

    int n;
    vector<vector<int>> adj;
    vector<Edge> E;
    vector<int> Q, ptr, dis, tree_parent, depth;
    vector<long long> tree_weight;

    GomoryHu(int n) : n(n), adj(n), Q(n), ptr(n), dis(n) {}

    /// One edge pair with capacity cap both ways: pushing f along one side frees f on the other
    void add_edge(int u, int v, long long cap){
        adj[u].push_back(E.size());
        E.push_back({v, cap, 0});

        adj[v].push_back(E.size());
        E.push_back({u, cap, 0});
    }

    void build(){
        tree_parent.assign(n, 0);
        tree_weight.assign(n, 0);

        for (int s = 1; s < n; s++){
            int t = tree_parent[s];
            for (auto& e : E) e.flow = 0;
            tree_weight[s] = maxflow(s, t);

            for (int v = 0; v < n; v++){
                if (v != s && tree_parent[v] == t && dis[v] != -1) tree_parent[v] = s;
            }

            /// Gusfield's swap: without it the tree is only flow-equivalent and subtrees are not min cuts
            if (dis[tree_parent[t]] != -1){
                tree_parent[s] = tree_parent[t], tree_parent[t] = s;
                swap(tree_weight[s], tree_weight[t]);
            }
        }

        vector<vector<int>> children(n);
        for (int v = 1; v < n; v++) children[tree_parent[v]].push_back(v);
        depth.assign(n, 0);
        vector<int> order = {0};
        for (int i = 0; i < (int)order.size(); i++){
            for (int c : children[order[i]]) depth[c] = depth[order[i]] + 1, order.push_back(c);
        }
    }

    long long min_cut(int u, int v) const{
        assert(u != v && (int)depth.size() == n);
        long long res = LLONG_MAX;
        while (u != v){
            if (depth[u] < depth[v]) swap(u, v);
            res = min(res, tree_weight[u]);
            u = tree_parent[u];
        }

        return res;
    }

    /// On return dis[v] != -1 exactly for the source side of a minimum cut, as the last bfs failed
    long long maxflow(int src, int sink){
        long long flow = 0;

        while (bfs(src, sink)){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, sink, LLONG_MAX)) flow += f;
        }

        return flow;
    }

    bool bfs(int src, int sink){
        int f = 0, l = 0;
        fill(dis.begin(), dis.end(), -1);

        dis[src] = 0, Q[l++] = src;
        while (f < l && dis[sink] == -1){
            int u = Q[f++];
            for (auto id : adj[u]){
                if (dis[E[id].v] == -1 && E[id].flow < E[id].cap){
                    Q[l++] = E[id].v;
                    dis[E[id].v] = dis[u] + 1;
                }
            }
        }

        return dis[sink] != -1;
    }

    long long dfs(int u, int sink, long long flow){
        if (u == sink || !flow) return flow;

        int len = adj[u].size();
        while (ptr[u] < len){
            int id = adj[u][ptr[u]];
            if (dis[E[id].v] == dis[u] + 1){
                long long f = dfs(E[id].v, sink, min(flow, E[id].cap - E[id].flow));
                if (f){
                    E[id].flow += f, E[id ^ 1].flow -= f;
                    return f;
                }
            }
            ptr[u]++;
        }

        return 0;
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

#define MAX 210

int n, cap[MAX][MAX];

int main(){
    int t, T, i, j;
    char buf[64];

    while (read(t)){
        T = 0;
        while (t--){
            read(n);
            for (i = 0; i < n; i++){
                for (j = 0; j < n; j++) read(cap[i][j]);
            }

            GomoryHu g(n);
            for (i = 0; i < n; i++){
                for (j = i + 1; j < n; j++){
                    if (cap[i][j]) g.add_edge(i, j, cap[i][j]);
                }
            }
            if (n) g.build();

            snprintf(buf, sizeof(buf), "Case #%d:", ++T);
            write(string(buf));
            for (i = 0; i < n; i++){
                vector<long long> row(n, 0);
                for (j = 0; j < n; j++){
                    if (i != j) row[j] = g.min_cut(i, j);
                }
                write(row);
            }
        }
    }

    return 0;
}
