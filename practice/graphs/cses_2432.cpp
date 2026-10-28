/***
 *
 * Maximum flow with Dinic's algorithm
 * Max flow and min cut on nodes 0 to n - 1, with optional node capacities or a dense capacity matrix
 *
 * Complexity: O(n^2 m) for m edges, O(m sqrt(n)) on unit-capacity bipartite matching
 *
 * FlowGraph g(n, src, sink): empty graph with edge capacities, src != sink
 * g.add_directed_edge(u, v, cap): adds u -> v and returns its id into g.E
 * g.add_undirected_edge(u, v, cap): adds u -> v and v -> u, both with capacity cap
 * g.maxflow(): pushes as much flow as possible and returns the flow this call added, the max flow on the first call
 * g.cut_edges(): ids into E of the edges leaving the source side, their capacities sum to the max flow
 * g.in_source_side(u): after maxflow(), true if src reaches u in the residual graph, the smallest min cut side
 * g.E[id].flow: flow on edge id after maxflow()
 *
 * FlowGraphWithNodeCap g(n, src, sink, node_capacity): the edge and maxflow() calls of FlowGraph, on g.flowgraph
 * DenseFlowGraph g(n, src, sink): an n x n capacity matrix, FlowGraph's calls without cut_edges() or per-edge flows
 *
 * Notes:
 *   - On bipartite matching hopcroft_karp.cpp is up to 2x faster on random graphs, on par on Library Checker's slowest
 *   - maxflow() keeps flat per-node copies of the edges, about 20 extra bytes per edge
 *   - maxflow() is the only entry point, bfs() and blocking_flow() are UB before it builds the flat copies
 *   - Before maxflow() the source side is empty and cut_edges() returns nothing
 *   - Node capacities apply to src and sink too, give them a huge capacity to leave them uncapped
 *   - FlowGraphWithNodeCap edge ids index g.flowgraph.E, where node x is 2x (in) and 2x + 1 (out)
 *   - Its min cut comes from g.flowgraph, a cut edge 2x -> 2x + 1 is a cut vertex x
 *   - Use DenseFlowGraph with half or more of all n^2 edges: 2x to 7x faster for n = 600 to 2000, 32 MB at n = 2000
 *   - Max closure, items with weights w[i] where picking i forces picking j:
 *       - src -> i with capacity w[i] if w[i] > 0, i -> sink with capacity -w[i] if w[i] < 0, i -> j with LLONG_MAX
 *       - The best total weight is the sum of positive w[i] - maxflow(), the picked items are those in_source_side(i)
 *       - The sum of positive w[i] must be below LLONG_MAX, the LLONG_MAX edges never overflow, in DenseFlowGraph too
 *   - Max density subgraph, max |E(S)| / |S| over nonempty S: Dinkelbach on the closure above, n + m + 2 nodes
 *       - Start with the density a / b = m / n of the whole graph
 *       - Closure items: each edge with weight b, forcing both endpoints, each vertex with weight -a
 *       - If the best closure is positive, set a / b to |E(S)| / |S| of its vertices S and repeat, else stop
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FlowGraph{
    struct Edge{
        int u, v;
        long long cap, flow;
    };

    int n, src, sink;
    vector<vector<int>> adj;
    vector<Edge> E;
    vector<int> Q, ptr, dis, path, start, to, rev, edge_id;
    vector<long long> res;

    FlowGraph(int n, int src, int sink): n(n), src(src), sink(sink), adj(n), Q(n), ptr(n), dis(n, -1) {}

    int add_directed_edge(int u, int v, long long cap){
        adj[u].push_back(E.size()), E.push_back({u, v, cap, 0});
        adj[v].push_back(E.size()), E.push_back({v, u, 0, 0});
        return E.size() - 2;
    }

    void add_undirected_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(){
        assert(src != sink);
        long long flow = 0;
        build_residual();

        while (bfs()){
            copy(start.begin(), start.end() - 1, ptr.begin());
            flow += blocking_flow();
        }

        for (int i = 0; i < (int)res.size(); i++) E[edge_id[i]].flow = E[edge_id[i]].cap - res[i];
        return flow;
    }

    vector<int> cut_edges() const{
        vector<int> ids;
        for (int id = 0; id < (int)E.size(); id++){
            if (E[id].cap > 0 && in_source_side(E[id].u) && !in_source_side(E[id].v)) ids.push_back(id);
        }

        return ids;
    }

    /// Reads dis from the final bfs of maxflow(), which failed to reach sink and so visited every node src reaches
    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 with residual res: reading E through adj missed the cache
    void build_residual(){
        int m = E.size();
        vector<int> pos(m);
        start.assign(n + 1, m), to.resize(m), rev.resize(m), edge_id.resize(m), res.resize(m);

        for (int u = 0, i = 0; u < n; u++){
            start[u] = i;
            for (int id : adj[u]) pos[id] = i, edge_id[i] = id, to[i] = E[id].v, res[i] = E[id].cap - E[id].flow, i++;
        }
        for (int i = 0; i < m; i++) rev[i] = pos[edge_id[i] ^ 1];
    }

    bool bfs(){
        int f = 0, l = 0;
        fill(dis.begin(), dis.end(), -1);

        dis[src] = 0, Q[l++] = src;
        while (f < l && dis[sink] == -1){
            int u = Q[f++];
            for (int i = start[u]; i < start[u + 1]; i++){
                if (res[i] > 0 && dis[to[i]] == -1) dis[to[i]] = dis[u] + 1, Q[l++] = to[i];
            }
        }

        return dis[sink] != -1;
    }

    /// Iterative walk from sink enters only nodes src reaches (2.5x-3x faster), path[j] is the arc taken j nodes back
    long long blocking_flow(){
        long long flow = 0;
        int u = sink;
        path.clear();

        while (true){
            if (u == src){
                long long f = LLONG_MAX;
                for (int i : path) f = min(f, res[rev[i]]);

                int k = path.size();
                for (int j = k - 1; j >= 0; j--){
                    res[rev[path[j]]] -= f, res[path[j]] += f;
                    if (res[rev[path[j]]] == 0) k = j;
                }
                flow += f, u = to[rev[path[k]]];
                path.resize(k);
                continue;
            }

            while (ptr[u] < start[u + 1] && (res[rev[ptr[u]]] == 0 || dis[to[ptr[u]]] != dis[u] - 1)) ptr[u]++;

            if (ptr[u] < start[u + 1]){
                path.push_back(ptr[u]);
                u = to[ptr[u]];
            }
            else{
                if (u == sink) return flow;
                u = to[rev[path.back()]];
                path.pop_back();
                ptr[u]++;
            }
        }
    }
};

struct FlowGraphWithNodeCap{
    FlowGraph flowgraph;

    FlowGraphWithNodeCap(int n, int src, int sink, vector<long long> node_capacity) : flowgraph(2 * n, 2 * src, 2 * sink + 1){
        assert(src != sink);
        for (int i = 0; i < n; i++){
            flowgraph.add_directed_edge(2 * i, 2 * i + 1, node_capacity[i]);
        }
    }

    int add_directed_edge(int u, int v, long long cap){
        return flowgraph.add_directed_edge(2 * u + 1, 2 * v, cap);
    }

    void add_undirected_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(){
        return flowgraph.maxflow();
    }
};

/// cap[u][v] is the residual capacity of u -> v, no flow is stored
struct DenseFlowGraph{
    int n, src, sink;
    vector<vector<long long>> cap;
    vector<int> dis, ptr;

    DenseFlowGraph(int n, int src, int sink) : n(n), src(src), sink(sink), cap(n, vector<long long>(n, 0)), dis(n, -1), ptr(n) {}

    void add_directed_edge(int u, int v, long long c){
        saturating_add(cap[u][v], c);
    }

    void add_undirected_edge(int u, int v, long long c){
        add_directed_edge(u, v, c);
        add_directed_edge(v, u, c);
    }

    long long maxflow(){
        assert(src != sink);
        long long flow = 0;

        while (bfs()){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, LLONG_MAX)) flow += f;
        }

        return flow;
    }

    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    /// One cell holds u -> v plus the flow on v -> u, so it saturates at LLONG_MAX, which still means infinite
    static void saturating_add(long long& cell, long long c){
        cell = cell > LLONG_MAX - c ? LLONG_MAX : cell + c;
    }

    bool bfs(){
        fill(dis.begin(), dis.end(), -1);
        vector<int> queue = {src};
        dis[src] = 0;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int v = 0; v < n; v++){
                if (dis[v] == -1 && cap[u][v] > 0) dis[v] = dis[u] + 1, queue.push_back(v);
            }
        }

        return dis[sink] != -1;
    }

    long long dfs(int u, long long flow){
        if (u == sink || !flow) return flow;
        for (int& v = ptr[u]; v < n; v++){
            if (dis[v] != dis[u] + 1 || cap[u][v] <= 0) continue;
            long long f = dfs(v, min(flow, cap[u][v]));
            if (f){
                cap[u][v] -= f, saturating_add(cap[v][u], f);
                return f;
            }
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

#define MAX 55

int id[MAX][MAX];

int main(){
    int n, i, j, x, sa = 0, sb = 0;

    read(n);
    FlowGraph g(2 * n + 2, 0, 2 * n + 1);
    for (i = 1; i <= n; i++) read(x), sa += x, g.add_directed_edge(0, i, x);
    for (i = 1; i <= n; i++) read(x), sb += x, g.add_directed_edge(n + i, 2 * n + 1, x);
    for (i = 1; i <= n; i++){
        for (j = 1; j <= n; j++) id[i][j] = g.add_directed_edge(i, n + j, 1);
    }

    if (sa != sb || g.maxflow() != sa){
        write(-1);
        return 0;
    }

    for (i = 1; i <= n; i++){
        string row(n, '.');
        for (j = 1; j <= n; j++){
            if (g.E[id[i][j]].flow > 0) row[j - 1] = 'X';
        }
        write(row);
    }

    return 0;
}
