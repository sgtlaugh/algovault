/***
 *
 * Maximum flow with highest-label push-relabel (HLPP)
 * Max flow and min cut on nodes 0 to n - 1, a drop-in for maxflow.cpp's FlowGraph
 *
 * Complexity: O(n^2 sqrt(m)) for m edges
 *
 * HLPP<T = long long> g(n, src, sink): empty graph with edge capacities of type T, src != sink
 * g.add_directed_edge(u, v, cap): adds u -> v and returns its id into g.E
 * g.add_undirected_edge(u, v, cap): adds u -> v and v -> u, both with capacity cap
 * g.maxflow(): pushes as much flow as possible and returns the flow this call added
 * g.cut_edges(): ids into E of the edges leaving the source side, their capacities sum to the max flow
 * g.in_source_side(u): after maxflow(), true if src reaches u in the residual graph, the smallest min cut side
 * g.E[id].flow: flow on edge id after maxflow()
 *
 * Notes:
 *   - Capacities of edges touching src must sum to at most the max of T, as maxflow() pushes them at once, asserted
 *   - Unlike FlowGraph, maxflow() takes no limit
 *   - Each maxflow() call continues from the flow already pushed, so it can be called again after adding edges
 *   - Faster than FlowGraph where Dinic needs many phases: 3x on a 300 x 300 grid, 7x on a 100-layer graph
 *   - On par with FlowGraph on random sparse and dense graphs, up to 2x slower on bipartite perfect matchings
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T = long long>
struct HLPP{
    struct Edge{
        int u, v;
        T cap, flow;
    };

    int n, src, sink;
    int highest, max_level, work;
    vector<vector<int>> adj;
    vector<Edge> E;
    vector<int> start, to, rev, edge_id;
    vector<int> height, cur, active, next_active, prev_level, next_level, order;
    vector<T> res, excess;
    vector<bool> side;

    HLPP(int n, int src, int sink): n(n), src(src), sink(sink), adj(n), side(n) {}

    int add_directed_edge(int u, int v, T cap){
        adj[u].push_back(E.size()), E.push_back({u, v, cap, 0});
        adj[v].push_back(E.size()), E.push_back({v, u, 0, 0});
        return E.size() - 2;
    }

    void add_undirected_edge(int u, int v, T cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    T maxflow(){
        assert(src != sink);
        build_residual();
        excess.assign(n, 0);

        T emitted = 0;
        for (int i = start[src]; i < start[src + 1]; i++){
            bool overflow = __builtin_add_overflow(emitted, res[i], &emitted);
            assert(!overflow);
            excess[to[i]] += res[i], res[rev[i]] += res[i], res[i] = 0;
        }

        /// The first pass leaves a max preflow, whose leftover excess the second pass returns to src
        push_relabel(sink, src);
        push_relabel(src, sink);

        for (int i = 0; i < (int)res.size(); i++) E[edge_id[i]].flow = E[edge_id[i]].cap - res[i];
        find_source_side();
        return excess[sink];
    }

    vector<int> cut_edges() const{
        vector<int> ids;
        for (int id = 0; id < (int)E.size(); id++){
            if (E[id].cap > 0 && in_source_side(E[id].u) && !in_source_side(E[id].v)) ids.push_back(id);
        }

        return ids;
    }

    bool in_source_side(int u) const{
        return side[u];
    }

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 with residual res, as in FlowGraph
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

    /// Moves every excess outside {target, blocked} into target; blocked keeps height n so nothing enters it
    void push_relabel(int target, int blocked){
        global_relabel(target, blocked);

        while (highest >= 0){
            int u = active[highest];
            if (u == -1){
                highest--;
                continue;
            }

            active[highest] = next_active[u];
            discharge(u, target);
            /// Global relabel once relabels scan 4n + m arcs: 1n to 8n plus m measure alike, 6n alone is 3x slower
            if (work > 4 * n + (int)res.size()) global_relabel(target, blocked);
        }
    }

    /// height[u] is the BFS distance from u to target in the residual graph, n when target is unreachable
    void global_relabel(int target, int blocked){
        height.assign(n, n), cur.resize(n), active.assign(n, -1), next_active.resize(n);
        prev_level.resize(2 * n), next_level.resize(2 * n);
        for (int k = n; k < 2 * n; k++) prev_level[k] = next_level[k] = k;

        height[target] = 0, order.assign(1, target);
        for (int k = 0; k < (int)order.size(); k++){
            int v = order[k];
            for (int i = start[v]; i < start[v + 1]; i++){
                int u = to[i];
                if (height[u] == n && u != blocked && res[rev[i]] > 0) height[u] = height[v] + 1, order.push_back(u);
            }
        }

        highest = -1, work = 0;
        for (int u : order){
            cur[u] = start[u], max_level = height[u];
            link_level(u);
            if (u != target && excess[u] > 0) activate(u), highest = height[u];
        }
    }

    /// Node n + k heads the circular list of the nodes at height k, which the gap heuristic empties in one sweep
    void link_level(int u){
        int head = n + height[u];
        prev_level[u] = head, next_level[u] = next_level[head];
        prev_level[next_level[head]] = u, next_level[head] = u;
    }

    void activate(int u){
        next_active[u] = active[height[u]], active[height[u]] = u;
    }

    void discharge(int u, int target){
        while (excess[u] > 0){
            if (cur[u] == start[u + 1]){
                relabel(u);
                if (height[u] == n) return;
                highest = height[u];
                continue;
            }

            int i = cur[u], v = to[i];
            if (res[i] > 0 && height[u] == height[v] + 1){
                T pushed = min(excess[u], res[i]);
                if (v != target && excess[v] == 0) activate(v);
                excess[u] -= pushed, excess[v] += pushed, res[i] -= pushed, res[rev[i]] += pushed;
            }
            else cur[u]++;
        }
    }

    /// Gap heuristic: when u leaves the last node at its height, every node above it can no longer reach target
    void relabel(int u){
        int old = height[u];
        prev_level[next_level[u]] = prev_level[u], next_level[prev_level[u]] = next_level[u];

        height[u] = n;
        for (int i = start[u]; i < start[u + 1]; i++){
            if (res[i] > 0 && height[to[i]] + 1 < height[u]) height[u] = height[to[i]] + 1, cur[u] = i;
        }
        work += start[u + 1] - start[u] + 12;

        if (next_level[n + old] == n + old){
            for (int k = old + 1; k <= max_level; k++){
                for (int v = next_level[n + k]; v != n + k; v = next_level[v]) height[v] = n;
                prev_level[n + k] = next_level[n + k] = n + k;
            }
            height[u] = n, max_level = old - 1;
        }
        else if (height[u] < n) link_level(u), max_level = max(max_level, height[u]);
    }

    void find_source_side(){
        fill(side.begin(), side.end(), false);
        order.assign(1, src);
        side[src] = true;

        for (int k = 0; k < (int)order.size(); k++){
            int u = order[k];
            for (int i = start[u]; i < start[u + 1]; i++){
                if (res[i] > 0 && !side[to[i]]) side[to[i]] = true, order.push_back(to[i]);
            }
        }
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

#define MAXC 300010

const int INF = 1000000000;

char buf[64];
int owner[MAXC];

vector<int> read_bid(){
    string line;
    vector<int> res;

    while (res.empty() && read_line(line)){
        int x = 0;
        bool digit = false;
        for (char c : line){
            if (c >= '0' && c <= '9') x = (x * 10) + (c - '0'), digit = true;
            else if (digit) res.push_back(x), x = 0, digit = false;
        }
        if (digit) res.push_back(x);
    }

    return res;
}

int main(){
    int t, tc, n, m, i, j, total;

    t = read_bid()[0];
    for (tc = 1; tc <= t; tc++){
        n = read_bid()[0];
        vector<vector<int>> telecom(n);
        for (i = 0; i < n; i++) telecom[i] = read_bid();

        m = read_bid()[0];
        vector<vector<int>> mobile(m);
        for (i = 0; i < m; i++) mobile[i] = read_bid();

        total = 0;
        HLPP<int> g(n + m + 2, 0, n + m + 1);
        for (i = 0; i < n; i++){
            g.add_directed_edge(0, i + 1, telecom[i][0]), total += telecom[i][0];
            for (j = 1; j < (int)telecom[i].size(); j++) owner[telecom[i][j]] = i + 1;
        }
        for (i = 0; i < m; i++){
            g.add_directed_edge(n + 1 + i, n + m + 1, mobile[i][0]), total += mobile[i][0];
            for (j = 1; j < (int)mobile[i].size(); j++){
                if (owner[mobile[i][j]]) g.add_directed_edge(owner[mobile[i][j]], n + 1 + i, INF);
            }
        }
        for (i = 0; i < n; i++){
            for (j = 1; j < (int)telecom[i].size(); j++) owner[telecom[i][j]] = 0;
        }

        if (tc > 1) write("");
        snprintf(buf, sizeof(buf), "Case %d:", tc);
        write(string(buf));
        write(total - g.maxflow());
    }

    return 0;
}
