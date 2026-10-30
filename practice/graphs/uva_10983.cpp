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
 * g.maxflow(limit = LLONG_MAX): pushes as much flow as possible up to limit and returns the flow this call added
 * g.cut_edges(): ids into E of the edges leaving the source side, their capacities sum to the max flow
 * g.in_source_side(u): after maxflow(), true if src reaches u in the residual graph, the smallest min cut side
 * g.E[id].flow: flow on edge id after maxflow()
 *
 * FlowGraphWithNodeCap g(n, src, sink, node_capacity): the edge and maxflow() calls of FlowGraph, on g.flowgraph
 * DenseFlowGraph g(n, src, sink): an n x n capacity matrix, FlowGraph's calls without cut_edges() or per-edge flows
 *
 * Notes:
 *   - On bipartite matching hopcroft_karp.cpp is up to 1.8x faster on random graphs, on par on Library Checker's worst
 *   - Memory is about 100 bytes per add_directed_edge call: two entries of E and maxflow()'s flat per-node copies
 *   - maxflow() is the only entry point, bfs() and blocking_flow() are UB before it builds the flat copies
 *   - Before maxflow() the source side is empty and cut_edges() returns nothing
 *   - A call stopped by its limit below the max flow leaves in_source_side(sink) true and cut_edges() no min cut
 *   - Each maxflow() call continues from the flow already pushed, so an unlimited call after a limited one completes it
 *   - Node capacities apply to src and sink too, give them a huge capacity to leave them uncapped
 *   - FlowGraphWithNodeCap edge ids index g.flowgraph.E, where node x is 2x (in) and 2x + 1 (out)
 *   - Its min cut comes from g.flowgraph, a cut edge 2x -> 2x + 1 is a cut vertex x
 *   - Use DenseFlowGraph with half or more of all n^2 edges: 3x to 10x faster for n = 600 to 2000, 32 MB at n = 2000
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
    vector<Edge> E;
    vector<int> Q, ptr, dis, path, start, to, rev, edge_id;
    vector<long long> res;

    FlowGraph(int n, int src, int sink): n(n), src(src), sink(sink), Q(n), ptr(n), dis(n, -1) {}

    int add_directed_edge(int u, int v, long long cap){
        E.push_back({u, v, cap, 0});
        E.push_back({v, u, 0, 0});
        return E.size() - 2;
    }

    void add_undirected_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(long long limit = LLONG_MAX){
        assert(src != sink);
        long long flow = 0;
        build_residual();

        /// bfs() runs before the limit check, so in_source_side(sink) tells whether the flow is maximum
        while (bfs() && flow < limit){
            copy(start.begin(), start.end() - 1, ptr.begin());
            flow += blocking_flow(limit - flow);
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

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 in id order, counting sorted: adjacency lists cost 2x on sparse
    void build_residual(){
        int m = E.size();
        vector<int> pos(m);

        start.assign(n + 1, 0), to.resize(m), rev.resize(m), edge_id.resize(m), res.resize(m);

        for (auto& e : E) start[e.u + 1]++;
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        vector<int> fill_at(start.begin(), start.end() - 1);
        for (int id = 0; id < m; id++){
            int i = pos[id] = fill_at[E[id].u]++;
            edge_id[i] = id, to[i] = E[id].v, res[i] = E[id].cap - E[id].flow;
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
    long long blocking_flow(long long limit){
        long long flow = 0;
        int u = sink;

        path.clear();

        while (true){
            if (u == src){
                long long f = limit - flow;
                for (int i : path) f = min(f, res[rev[i]]);

                int k = path.size();
                for (int j = k - 1; j >= 0; j--){
                    res[rev[path[j]]] -= f, res[path[j]] += f;
                    if (res[rev[path[j]]] == 0) k = j;
                }
                flow += f;
                if (flow == limit) return flow;

                u = to[rev[path[k]]];
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

    long long maxflow(long long limit = LLONG_MAX){
        return flowgraph.maxflow(limit);
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

    long long maxflow(long long limit = LLONG_MAX){
        assert(src != sink);
        long long flow = 0;

        while (bfs() && flow < limit){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, limit - flow)) flow += f;
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
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
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

#define MAXM 1010

int n, d, m, total, U[MAXM], V[MAXM], C[MAXM], P[MAXM], E[MAXM], Z[35];

int node(int city, int day){
    return day * n + city;
}

bool possible(int cost){
    int i, src = n * (d + 1), sink = node(n - 1, d);

    FlowGraph g(src + 1, src, sink);
    for (i = 0; i < n; i++) g.add_directed_edge(src, node(i, 0), Z[i]);
    for (i = 0; i < (n * d); i++) g.add_directed_edge(i, i + n, total);
    for (i = 0; i < m; i++){
        if (P[i] <= cost) g.add_directed_edge(node(U[i], E[i]), node(V[i], E[i] + 1), C[i]);
    }

    return g.maxflow() == total;
}

int main(){
    int t, T, i, low, high, mid;
    char buf[64];

    while (read(t)){
        T = 0;
        while (t--){
            read(n, d, m);
            for (i = 0; i < m; i++){
                read(U[i], V[i], C[i], P[i], E[i]);
                U[i]--, V[i]--;
            }

            total = 0;
            for (i = 0; i < n; i++) read(Z[i]), total += Z[i];

            vector<int> costs = {0};
            for (i = 0; i < m; i++) costs.push_back(P[i]);
            sort(costs.begin(), costs.end());

            if (!possible(costs.back())) snprintf(buf, sizeof(buf), "Case #%d: Impossible", ++T);
            else{
                low = 0, high = costs.size() - 1;
                while (low < high){
                    mid = (low + high) >> 1;
                    if (possible(costs[mid])) high = mid;
                    else low = mid + 1;
                }
                snprintf(buf, sizeof(buf), "Case #%d: %d", ++T, costs[low]);
            }
            write(string(buf));
        }
    }

    return 0;
}
