/***
 *
 * Network Simplex
 * Cheapest flow meeting every node's supply or demand, on edges with lower and upper bounds and costs of either sign
 *
 * Complexity: O(n + m) per pivot, with no polynomial bound on the number of pivots
 *
 * NetworkSimplex g(n): empty graph on nodes 0 to n - 1, every supply 0
 * g.add_edge(u, v, cap, cost): adds u -> v carrying 0 to cap units, returns its id
 * g.add_edge(u, v, low, high, cost): adds u -> v carrying low to high units, low <= high of either sign, returns its id
 * g.add_supply(v, amount): node v must send amount more units out than in, a negative amount is a demand
 * g.solve(): the minimum total cost, nullopt when no flow meets every bound and supply
 * g.flow(id): flow on edge id after a feasible solve
 * g.potential[v]: optimal dual, u -> v has cost + potential[u] - potential[v] >= 0 if flow < high, <= 0 if flow > low
 *
 * Notes:
 *   - Requires __int128 for the total cost
 *   - Self loops, parallel edges and negative cycles are allowed
 *   - n max|cost| <= 1e18, and then |potential[v]| <= 2 n max|cost| + 1
 *   - sum |supply| + sum |low| + sum (high - low) must fit in long long
 *   - Several times faster than min_cost_circulation.cpp in practice, which instead guarantees a polynomial bound
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct NetworkSimplex{
    struct Edge{
        int from, to;
        long long low, cap, cost, flow;
    };

    static constexpr long long INF = numeric_limits<long long>::max();

    int n, m = 0;
    vector<Edge> edges;
    vector<long long> supply, potential;
    vector<int> parent, parent_edge, depth, first_child, next_sibling, prev_sibling, pending;
    int next_edge = 0;

    NetworkSimplex(int n) : n(n), supply(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        return add_edge(u, v, 0, cap, cost);
    }

    /// cap and flow are stored shifted down by low, so the simplex only sees 0 <= flow <= cap
    int add_edge(int u, int v, long long low, long long high, long long cost){
        assert(0 <= u && u < n && 0 <= v && v < n && low <= high);
        edges.resize(m);
        edges.push_back({u, v, low, high - low, cost, 0});
        return m++;
    }

    void add_supply(int v, long long amount){
        assert(0 <= v && v < n);
        supply[v] += amount;
    }

    optional<__int128> solve(){
        int root = n;
        edges.resize(m);
        vector<long long> b = supply;
        long long max_cost = 1;
        for (auto& e : edges){
            e.flow = 0;
            b[e.from] -= e.low, b[e.to] += e.low;
            max_cost = max(max_cost, abs(e.cost));
        }

        /// A cycle through the root uses two root arcs, so big > n max_cost / 2 empties them when a flow is feasible
        long long big = n * max_cost + 1;
        potential.assign(n + 1, 0);
        parent.assign(n + 1, root), parent_edge.assign(n + 1, -1), depth.assign(n + 1, 1);
        first_child.assign(n + 1, -1), next_sibling.assign(n + 1, -1), prev_sibling.assign(n + 1, -1);
        parent[root] = -1, depth[root] = 0;

        /// Root arcs point so a zero flow one leads to the root: the tree starts strongly feasible, preventing cycling
        for (int v = 0; v < n; v++){
            if (b[v] >= 0) edges.push_back({v, root, 0, INF, big, b[v]}), potential[v] = -big;
            else edges.push_back({root, v, 0, INF, big, -b[v]}), potential[v] = big;
            parent_edge[v] = m + v;
            attach(v, root);
        }

        next_edge = 0;
        int block = max(10, (int)sqrt(edges.size()));
        while (true){
            int enter = entering_edge(block);
            if (enter == -1) break;
            pivot(enter);
        }

        /// Also catches a nonzero total supply: empty root arcs would mean conservation at every node
        for (int v = 0; v < n; v++){
            if (edges[m + v].flow > 0) return nullopt;
        }

        __int128 total = 0;
        for (int id = 0; id < m; id++) total += (__int128)edges[id].cost * flow(id);
        potential.resize(n);
        return total;
    }

    long long flow(int id) const{
        return edges[id].low + edges[id].flow;
    }

    long long reduced_cost(const Edge& e) const{
        return e.cost + potential[e.from] - potential[e.to];
    }

    /// Block search: scans block edges at a time from where the last search stopped, takes the worst violation found
    int entering_edge(int block){
        int total = edges.size(), enter = -1, in_block = 0;
        long long worst = 0;

        for (int scanned = 0; scanned < total; scanned++){
            const Edge& e = edges[next_edge];
            long long rc = reduced_cost(e), violation = 0;
            if (rc < 0 && e.flow < e.cap) violation = -rc;
            if (rc > 0 && e.flow > 0) violation = rc;
            if (violation > worst) worst = violation, enter = next_edge;

            next_edge = next_edge + 1 == total ? 0 : next_edge + 1;
            if (++in_block == block){
                if (enter != -1) break;
                in_block = 0;
            }
        }

        return enter;
    }

    /// Pushes flow around the cycle of enter and the tree, the last blocking arc after the apex leaves the tree
    void pivot(int enter){
        Edge& e = edges[enter];
        bool forward = reduced_cost(e) < 0;
        int p = forward ? e.from : e.to, q = forward ? e.to : e.from;

        int apex_p = p, apex_q = q;
        while (apex_p != apex_q){
            if (depth[apex_p] >= depth[apex_q]) apex_p = parent[apex_p];
            else apex_q = parent[apex_q];
        }
        int apex = apex_p;

        /// The cycle runs apex -> p, enter, q -> apex: strict on the p side and <= after keeps the last minimum
        long long delta = INF;
        int leave = -1;
        bool leave_on_p = false;
        for (int x = p; x != apex; x = parent[x]){
            long long r = residual(x, true);
            if (r < delta) delta = r, leave = x, leave_on_p = true;
        }
        long long r_enter = forward ? e.cap - e.flow : e.flow;
        if (r_enter <= delta) delta = r_enter, leave = -1;
        for (int x = q; x != apex; x = parent[x]){
            long long r = residual(x, false);
            if (r <= delta) delta = r, leave = x, leave_on_p = false;
        }

        if (delta > 0){
            e.flow += forward ? delta : -delta;
            for (int x = p; x != apex; x = parent[x]){
                Edge& t = edges[parent_edge[x]];
                t.flow += t.to == x ? delta : -delta;
            }
            for (int x = q; x != apex; x = parent[x]){
                Edge& t = edges[parent_edge[x]];
                t.flow += t.from == x ? delta : -delta;
            }
        }
        if (leave == -1) return;

        /// The path from the near end of enter up to leave flips direction, then hangs from the far end through enter
        int x = leave_on_p ? p : q, new_parent = leave_on_p ? q : p, new_edge = enter;
        while (true){
            int old_parent = parent[x], old_edge = parent_edge[x];
            detach(x);
            attach(x, new_parent);
            parent_edge[x] = new_edge;
            if (x == leave) break;
            new_parent = x, new_edge = old_edge, x = old_parent;
        }
        update_subtree(leave_on_p ? p : q);
    }

    /// Room to push along the tree arc between x and its parent, downward toward x or upward away from it
    long long residual(int x, bool down) const{
        const Edge& t = edges[parent_edge[x]];
        return (t.to == x) == down ? t.cap - t.flow : t.flow;
    }

    void detach(int x){
        if (prev_sibling[x] != -1) next_sibling[prev_sibling[x]] = next_sibling[x];
        else first_child[parent[x]] = next_sibling[x];
        if (next_sibling[x] != -1) prev_sibling[next_sibling[x]] = prev_sibling[x];
    }

    void attach(int x, int new_parent){
        parent[x] = new_parent, prev_sibling[x] = -1, next_sibling[x] = first_child[new_parent];
        if (first_child[new_parent] != -1) prev_sibling[first_child[new_parent]] = x;
        first_child[new_parent] = x;
    }

    /// Tree arcs have reduced cost 0, so each potential follows from its parent's
    void update_subtree(int top){
        pending.assign(1, top);
        while (!pending.empty()){
            int x = pending.back();
            pending.pop_back();
            const Edge& t = edges[parent_edge[x]];
            potential[x] = potential[parent[x]] + (t.to == x ? t.cost : -t.cost);
            depth[x] = depth[parent[x]] + 1;
            for (int c = first_child[x]; c != -1; c = next_sibling[c]) pending.push_back(c);
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

int main(){
    int n, m, i, u, v, low, high;

    read(n, m);
    NetworkSimplex g(n);
    for (i = 0; i < m; i++){
        read(u, v, low, high);
        g.add_edge(u - 1, v - 1, low, high, 0);
    }

    if (!g.solve()){
        write("NO");
        return 0;
    }

    write("YES");
    for (i = 0; i < m; i++) write(g.flow(i));
    return 0;
}
