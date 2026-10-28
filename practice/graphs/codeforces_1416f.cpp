/***
 *
 * Flow with lower bounds (demands) on Dinic's algorithm
 * Every edge carries a flow f with low <= f <= high: feasible circulation, feasible s-t flow, max and min s-t flow
 *
 * Complexity: one Dinic run per call on n + 2 nodes and m + n + 2 edges, O(n^2 m)
 *
 * LowerBoundFlow g(n): empty graph on nodes 0 to n - 1
 * g.add_edge(u, v, low, high): directed edge u -> v, 0 <= low <= high, returns its index
 * g.circulation(): true if some flow satisfies every bound with conservation at every node
 * g.feasible_flow(s, t): value of some s-t flow satisfying every bound, conservation everywhere but s and t
 * g.max_flow(s, t): the largest value of an s-t flow satisfying every bound
 * g.min_flow(s, t): the smallest value of an s-t flow satisfying every bound
 * g.flow(id): flow on edge id after a call that succeeded
 *
 * Notes:
 *   - The three s-t calls return LowerBoundFlow::INFEASIBLE when no flow satisfies the bounds
 *   - The value is the net flow out of s and can be negative when edges into s have lower bounds
 *   - Calls never change the stored edges, so one graph answers any number of queries
 *   - The sum of all upper bounds must fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LowerBoundFlow{
    static constexpr long long INFEASIBLE = numeric_limits<long long>::min();

    struct Bound{
        int u, v;
        long long low, high;
    };

    struct Arc{
        int to;
        long long cap;
    };

    int n;
    vector<Bound> bounds;
    vector<Arc> arcs;
    vector<vector<int>> adj;
    vector<int> dis, ptr;

    LowerBoundFlow(int n) : n(n), adj(n + 2), dis(n + 2), ptr(n + 2) {}

    int add_edge(int u, int v, long long low, long long high){
        assert(0 <= u && u < n && 0 <= v && v < n && 0 <= low && low <= high);
        bounds.push_back({u, v, low, high});
        return bounds.size() - 1;
    }

    bool circulation(){
        return build();
    }

    long long feasible_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;
        return detach_terminals();
    }

    long long max_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;

        long long value = detach_terminals();
        return value + dinic(s, t);
    }

    long long min_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;

        long long value = detach_terminals();
        return value - dinic(t, s);
    }

    /// Arc 2 * id is edge id's residual arc and arc 2 * id + 1 its reverse, whose capacity is the flow above low
    long long flow(int id) const{
        assert(0 <= id && 2 * id + 1 < (int)arcs.size());
        return bounds[id].low + arcs[2 * id + 1].cap;
    }

    /// Feeds lower bounds from S = n to T = n + 1, with the s-t arcs last so detach_terminals finds them at the back
    bool build(int s = -1, int t = -1){
        int super_src = n, super_sink = n + 1;
        arcs.clear();
        for (auto& list : adj) list.clear();

        vector<long long> excess(n, 0);
        long long total_high = 0;
        for (auto& b : bounds){
            add_arc(b.u, b.v, b.high - b.low);
            excess[b.v] += b.low, excess[b.u] -= b.low;
            total_high += b.high;
        }

        long long demand = 0;
        for (int x = 0; x < n; x++){
            if (excess[x] > 0) add_arc(super_src, x, excess[x]), demand += excess[x];
            if (excess[x] < 0) add_arc(x, super_sink, -excess[x]);
        }
        if (s != -1) add_arc(t, s, total_high), add_arc(s, t, total_high);

        return dinic(super_src, super_sink) == demand;
    }

    void add_arc(int u, int v, long long cap){
        adj[u].push_back(arcs.size()), arcs.push_back({v, cap});
        adj[v].push_back(arcs.size()), arcs.push_back({u, 0});
    }

    long long dinic(int src, int sink){
        long long flow = 0;

        while (bfs(src, sink)){
            fill(ptr.begin(), ptr.end(), 0);
            flow += blocking_flow(src, sink);
        }

        return flow;
    }

    bool bfs(int src, int sink){
        fill(dis.begin(), dis.end(), -1);
        vector<int> queue = {src};
        dis[src] = 0;

        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int id : adj[u]){
                if (dis[arcs[id].to] == -1 && arcs[id].cap > 0) dis[arcs[id].to] = dis[u] + 1, queue.push_back(arcs[id].to);
            }
        }

        return dis[sink] != -1;
    }

    /// Iterative: level graphs as deep as n would overflow an 8 MB stack with a recursive DFS
    long long blocking_flow(int src, int sink){
        long long total = 0;
        vector<int> path;
        int u = src;

        while (true){
            if (u == sink){
                long long f = LLONG_MAX;
                for (int id : path) f = min(f, arcs[id].cap);
                for (int id : path) arcs[id].cap -= f, arcs[id ^ 1].cap += f;
                total += f;

                int keep = 0;
                while (arcs[path[keep]].cap > 0) keep++;
                path.resize(keep);
                u = keep ? arcs[path.back()].to : src;
                continue;
            }

            int& i = ptr[u];
            while (i < (int)adj[u].size() && (arcs[adj[u][i]].cap <= 0 || dis[arcs[adj[u][i]].to] != dis[u] + 1)) i++;
            if (i < (int)adj[u].size()){
                path.push_back(adj[u][i]);
                u = arcs[adj[u][i]].to;
                continue;
            }

            if (path.empty()) break;
            u = arcs[path.back() ^ 1].to;
            path.pop_back();
            ptr[u]++;
        }

        return total;
    }

    /// Removes the t -> s and s -> t arcs, returning their value; S and T stay saturated so no path crosses them
    long long detach_terminals(){
        int m = arcs.size();
        long long value = arcs[m - 3].cap - arcs[m - 1].cap;
        for (int id = m - 4; id < m; id++) arcs[id].cap = 0;
        return value;
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

const int dx[] = {0, 0, 1, -1}, dy[] = {-1, 1, 0, 0};
const char dir[] = "LRDU";

int n, m;
vector<vector<int>> ar, cost;
vector<string> res;

bool inside(int x, int y){
    return x >= 0 && x < n && y >= 0 && y < m;
}

bool critical(int x, int y){
    for (int d = 0; d < 4; d++){
        int a = x + dx[d], b = y + dy[d];
        if (inside(a, b) && ar[a][b] < ar[x][y]) return false;
    }

    return true;
}

bool solve(){
    int i, j, d, a, b, s = n * m, t = s + 1;
    LowerBoundFlow g(n * m + 2);
    vector<pair<int, int>> ids;
    vector<int> eid;

    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            int need = critical(i, j);
            if ((i + j) & 1){
                g.add_edge(i * m + j, t, need, 1);
                continue;
            }

            g.add_edge(s, i * m + j, need, 1);
            for (d = 0; d < 4; d++){
                a = i + dx[d], b = j + dy[d];
                if (inside(a, b) && ar[a][b] == ar[i][j]){
                    eid.push_back(g.add_edge(i * m + j, a * m + b, 0, 1));
                    ids.push_back({i * m + j, d});
                }
            }
        }
    }

    if (g.feasible_flow(s, t) == LowerBoundFlow::INFEASIBLE) return false;

    cost.assign(n, vector<int>(m, 0));
    res.assign(n, string(m, '?'));
    for (i = 0; i < (int)eid.size(); i++){
        if (!g.flow(eid[i])) continue;

        int x = ids[i].first / m, y = ids[i].first % m, d = ids[i].second;
        a = x + dx[d], b = y + dy[d];
        cost[x][y] = 1, cost[a][b] = ar[x][y] - 1;
        res[x][y] = dir[d], res[a][b] = dir[d ^ 1];
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            if (res[i][j] != '?') continue;
            for (d = 0; d < 4; d++){
                a = i + dx[d], b = j + dy[d];
                if (inside(a, b) && ar[a][b] < ar[i][j]){
                    cost[i][j] = ar[i][j] - ar[a][b], res[i][j] = dir[d];
                    break;
                }
            }
        }
    }

    return true;
}

int main(){
    int T, i, j;

    read(T);
    while (T--){
        read(n, m);
        ar.assign(n, vector<int>(m));
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++) read(ar[i][j]);
        }

        if (!solve()){
            write("NO");
            continue;
        }

        write("YES");
        for (i = 0; i < n; i++) write(cost[i]);
        for (i = 0; i < n; i++){
            vector<char> row(res[i].begin(), res[i].end());
            write(row);
        }
    }

    return 0;
}
