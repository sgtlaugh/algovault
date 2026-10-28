/***
 *
 * Dijkstra's Algorithm and 0-1 BFS
 * Shortest paths from one or more sources with non-negative edge weights, plus shortest path reconstruction
 *
 * Complexity:
 *   - Dijkstra: O((n + m) log m) per run with a binary heap
 *   - DenseDijkstra: O(n^2) per run, O(n^2) memory for the weight matrix (32 MB at n = 2000)
 *   - ZeroOneBFS: O(n + m) per run with a deque
 *
 * Dijkstra, DenseDijkstra and ZeroOneBFS share one API on 0-based vertices, shown for Dijkstra:
 * Dijkstra g(n): an empty graph on n vertices
 * g.add_edge(u, v, w): adds the directed edge u -> v of weight w, pass directed = false for an undirected one
 * g.run(src): fills g.dist and g.parent with the shortest paths from src
 * g.run(src, target): stops as soon as target is settled, so only dist[target] and path(target) are final
 * g.run(sources): the same from the nearest of several sources, all at distance 0, and also takes a target
 * g.path(v): the vertices of one shortest path from a source to v, empty if v is unreachable
 * g.dist[v]: the shortest distance from a source to v, or DIJKSTRA_INF if v is unreachable
 *
 * Notes:
 *   - Weights are long long in [0, DIJKSTRA_INF), 0 or 1 for ZeroOneBFS
 *   - Every shortest distance must stay below DIJKSTRA_INF (~2.3e18)
 *   - Parallel edges and self loops are allowed, and add_edge and run can be called in any order
 *   - DenseDijkstra beats Dijkstra once m passes about n^2 / 2, or when most edges relax, like on a complete DAG
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long DIJKSTRA_INF = numeric_limits<long long>::max() / 4;

vector<int> trace_path(const vector<long long>& dist, const vector<int>& parent, int v){
    if (dist[v] == DIJKSTRA_INF) return {};

    vector<int> path;
    for (; v != -1; v = parent[v]) path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

struct Dijkstra{
    int n;
    vector<tuple<int, int, long long>> arcs;
    vector<int> start;
    vector<pair<int, long long>> adj;  /// u's arcs are adj[start[u] .. start[u + 1]), rebuilt by run after add_edge
    vector<long long> dist;
    vector<int> parent;

    Dijkstra(int n) : n(n), start(n + 1), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w >= 0);
        arcs.push_back({u, v, w});
        if (!directed) arcs.push_back({v, u, w});
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        if (adj.size() != arcs.size()) build_adjacency();
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;

        for (int src : sources) dist[src] = 0, heap.push({0, src});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != dist[u]) continue;
            if (u == target) break;
            for (auto it = adj.begin() + start[u], end = adj.begin() + start[u + 1]; it != end; it++){
                auto [v, w] = *it;
                if (du + w < dist[v]) dist[v] = du + w, parent[v] = u, heap.push({dist[v], v});
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }

    void build_adjacency(){
        fill(start.begin(), start.end(), 0);
        for (auto [u, v, w] : arcs) start[u + 1]++;
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        adj.resize(arcs.size());
        vector<int> next_slot(start.begin(), start.end() - 1);
        for (auto [u, v, w] : arcs) adj[next_slot[u]++] = {v, w};
    }
};

struct DenseDijkstra{
    int n;
    vector<long long> weight;  /// weight[u * n + v], the lightest u -> v edge or DIJKSTRA_INF
    vector<long long> dist;
    vector<int> parent;

    DenseDijkstra(int n) : n(n), weight((size_t)n * n, DIJKSTRA_INF), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w >= 0);
        long long& cell = weight[(size_t)u * n + v];
        cell = min(cell, w);
        if (!directed) add_edge(v, u, w, true);
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        vector<bool> done(n, false);

        for (int src : sources) dist[src] = 0;
        for (int it = 0; it < n; it++){
            int u = -1;
            for (int v = 0; v < n; v++){
                if (!done[v] && (u == -1 || dist[v] < dist[u])) u = v;
            }
            if (dist[u] == DIJKSTRA_INF || u == target) break;

            done[u] = true;
            /// A missing edge holds DIJKSTRA_INF, so dist[u] + INF >= INF >= dist[v] never relaxes and needs no branch
            const long long* row = &weight[(size_t)u * n];
            for (int v = 0; v < n; v++){
                if (dist[u] + row[v] < dist[v]) dist[v] = dist[u] + row[v], parent[v] = u;
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }
};

struct ZeroOneBFS{
    int n;
    vector<vector<pair<int, int>>> adj;
    vector<long long> dist;
    vector<int> parent;

    ZeroOneBFS(int n) : n(n), adj(n), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w == 0 || w == 1);
        adj[u].push_back({v, (int)w});
        if (!directed) adj[v].push_back({u, (int)w});
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        deque<int> dq;

        /// The deque only ever holds distances d and d + 1, so a vertex improves at most once after it is first reached
        for (int src : sources) dist[src] = 0, dq.push_back(src);
        while (!dq.empty()){
            int u = dq.front();
            dq.pop_front();
            if (u == target) break;
            for (auto [v, w] : adj[u]){
                if (dist[u] + w >= dist[v]) continue;
                dist[v] = dist[u] + w, parent[v] = u;
                if (w == 0) dq.push_front(v);
                else dq.push_back(v);
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
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

int main(){
    int T, n, m, t, s, g, h, i, a, b, d, x;
    long long gh = 0;

    read(T);
    while (T--){
        read(n, m, t);
        read(s, g, h);
        Dijkstra dij(n + 1);
        for (i = 0; i < m; i++){
            read(a, b, d);
            dij.add_edge(a, b, d, false);
            if ((a == g && b == h) || (a == h && b == g)) gh = d;
        }

        dij.run(s);
        vector<long long> ds = dij.dist;
        dij.run(g);
        vector<long long> dg = dij.dist;
        dij.run(h);
        vector<long long> dh = dij.dist;

        vector<int> res;
        for (i = 0; i < t; i++){
            read(x);
            if (ds[x] == DIJKSTRA_INF) continue;
            if (ds[g] + gh + dh[x] == ds[x] || ds[h] + gh + dg[x] == ds[x]) res.push_back(x);
        }

        sort(res.begin(), res.end());
        write(res);
    }

    return 0;
}
