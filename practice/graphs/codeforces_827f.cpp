/***
 *
 * Dijkstra's Algorithm and 0-1 BFS
 * Shortest paths from one or more sources with non-negative edge weights, plus shortest path reconstruction
 *
 * Complexity:
 *   - Dijkstra: O((n + m) log m) per run with a binary heap
 *   - DenseDijkstra: O(n^2) per run, O(n^2) memory for the weight matrix (32 MB at n = 2000)
 *   - ZeroOneBFS: O(n + m) per run with a deque
 *   - dijkstra_implicit: O((n + m) log m) for m relax calls, plus the cost of neighbors
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
 * dijkstra_implicit(n, sources, neighbors): the distances over states 0 .. n - 1, DIJKSTRA_INF where unreached
 * dijkstra_implicit(n, sources, neighbors, target): stops once target is settled, so only its distance is final
 * neighbors(u, du, DijkstraRelax& relax): calls relax(v, dv) for each successor v of u, with any dv >= du
 *
 * Notes:
 *   - Weights are long long in [0, DIJKSTRA_INF), 0 or 1 for ZeroOneBFS
 *   - Every shortest distance must stay below DIJKSTRA_INF (~2.3e18)
 *   - Parallel edges and self loops are allowed, and add_edge and run can be called in any order
 *   - DenseDijkstra beats Dijkstra once m passes about n^2 / 2, or when most edges relax, like on a complete DAG
 *   - dijkstra_implicit's dv may depend on du, like a toll on the load carried, but must never drop as du grows
 *   - dijkstra_implicit calls neighbors at most once per state, with its final du
 *   - A start distance other than 0 comes from an extra source state whose neighbors relax the real start to it
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

struct DijkstraRelax{
    vector<long long>& dist;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>>& heap;
    long long du;

    void operator()(int v, long long dv){
        assert(dv >= du);
        if (dv < dist[v]) dist[v] = dv, heap.push({dv, v});
    }
};

template<typename Neighbors>
vector<long long> dijkstra_implicit(int n, const vector<int>& sources, Neighbors neighbors, int target = -1){
    vector<long long> dist(n, DIJKSTRA_INF);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
    DijkstraRelax relax{dist, heap, 0};

    /// A repeated source is pushed once, or neighbors would run twice for it
    for (int src : sources) if (dist[src] != 0) dist[src] = 0, heap.push({0, src});
    while (!heap.empty()){
        auto [du, u] = heap.top();
        heap.pop();
        if (du != dist[u]) continue;
        if (u == target) break;
        relax.du = du;
        neighbors(u, du, relax);
    }

    return dist;
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
        for (; isdigit(c); c = read_char()){
            v = v * 10 - (c - '0');
        }
        x = neg ? v : -v;

        return true;
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return 0;

        for (; !isspace(c) && c != EOF; c = read_char()){
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

#define MAX 500010

int n, m, to[2 * MAX], lo[2 * MAX], hi[2 * MAX], ptr[MAX][2];
vector<int> out[MAX];

void neighbors(int u, long long t, DijkstraRelax& relax){
    int e = u >> 1, q = (t + 1) & 1, b = to[e];
    long long R = hi[e] - ((hi[e] - t - 1) & 1);

    while (ptr[b][q] < (int)out[b].size()){
        int f = out[b][ptr[b][q]];
        long long x = max(t + 1, (long long)lo[f]);
        if ((x & 1) != q) x++;
        if (x > R) break;

        if (x < hi[f]) relax(2 * f + q, x);
        ptr[b][q]++;
    }
}

int main(){
    int i, a, b, l, r;

    read(n, m);
    for (i = 0; i < m; i++){
        read(a, b, l, r);
        to[2 * i] = b, to[2 * i + 1] = a;
        lo[2 * i] = lo[2 * i + 1] = l;
        hi[2 * i] = hi[2 * i + 1] = r;
        out[a].push_back(2 * i);
        out[b].push_back(2 * i + 1);
    }

    if (n == 1){
        write(0);
        return 0;
    }

    vector<int> sources;
    for (i = 1; i <= n; i++) sort(out[i].begin(), out[i].end(), [](int x, int y){ return lo[x] < lo[y]; });
    for (int e: out[1]){
        if (lo[e] == 0) sources.push_back(2 * e);
    }

    vector<long long> dist = dijkstra_implicit(4 * m, sources, neighbors);

    long long res = DIJKSTRA_INF;
    for (i = 0; i < 4 * m; i++){
        if (to[i >> 1] == n) res = min(res, dist[i] + 1);
    }

    write(res >= DIJKSTRA_INF ? -1 : res);
    return 0;
}
