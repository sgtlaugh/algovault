/***
 *
 * Minimum Cost Circulation
 * Cheapest flow with zero net flow at every node, negative costs and negative cycles allowed
 *
 * Complexity: O(n^3 log(n C)) for C = max |cost|
 *
 * MinCostCirculation g(n): empty graph on nodes 0 to n - 1
 * g.add_edge(u, v, cap, cost): adds u -> v with cap >= 0 and a cost of either sign, returns its id
 * g.solve(): the minimum cost, 0 when no negative cycle exists
 * g.flow(id): flow on edge id in that optimal circulation
 *
 * Notes:
 *   - Self loops and parallel edges are allowed
 *   - n <= 1e4 and |cost| <= 1e9
 *   - The total cost sum(flow * cost) and the capacity sum into any node must fit in long long
 *   - Min cost max flow from s to t with negative cycles: add_edge(t, s, sum of caps out of s, -M), M = 1 + sum |cost|
 *       - Max flow is the flow on that edge, min cost is solve() + M * max flow, M must respect the cost limit above
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MinCostCirculation{
    struct Arc{
        int to, rev;
        long long cap, cost;
    };

    static constexpr long long ALPHA = 4;

    int n;
    vector<vector<Arc>> adj;
    vector<pair<int, int>> pos;
    vector<long long> excess, potential;
    vector<int> cur;

    MinCostCirculation(int n) : n(n), adj(n), excess(n), potential(n), cur(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        assert(0 <= u && u < n && 0 <= v && v < n && cap >= 0);
        pos.push_back({u, (int)adj[u].size()});

        /// Scaled by n + 1: a cycle has at most n arcs, so 1-optimality in scaled units leaves no negative cycle
        cost *= n + 1;
        adj[u].push_back({v, (int)adj[v].size() + (u == v), cap, cost});
        adj[v].push_back({u, (int)adj[u].size() - 1, 0, -cost});

        return pos.size() - 1;
    }

    long long solve(){
        long long eps = 0;
        for (auto& arcs : adj){
            for (auto& a : arcs) eps = max(eps, abs(a.cost));
        }

        fill(potential.begin(), potential.end(), 0);
        while (eps > 1){
            eps = max(1LL, eps / ALPHA);
            refine(eps);
        }

        long long total = 0;
        for (int id = 0; id < (int)pos.size(); id++){
            auto [u, i] = pos[id];
            total += flow(id) * (adj[u][i].cost / (n + 1));
        }
        return total;
    }

    long long flow(int id) const{
        auto [u, i] = pos[id];
        const Arc& a = adj[u][i];
        return adj[a.to][a.rev].cap;
    }

    /// Turns the circulation eps-optimal: saturate every negative reduced cost arc, then discharge the excesses FIFO
    void refine(long long eps){
        for (int u = 0; u < n; u++){
            for (auto& a : adj[u]){
                if (a.cap > 0 && reduced_cost(u, a) < 0) push(u, a, a.cap);
            }
        }

        queue<int> active;
        vector<bool> queued(n, false);
        for (int u = 0; u < n; u++){
            cur[u] = 0;
            if (excess[u] > 0) active.push(u), queued[u] = true;
        }

        while (!active.empty()){
            int u = active.front();
            active.pop(), queued[u] = false;
            while (excess[u] > 0){
                if (cur[u] == (int)adj[u].size()) relabel(u, eps);
                Arc& a = adj[u][cur[u]];
                if (a.cap == 0 || reduced_cost(u, a) >= 0){
                    cur[u]++;
                    continue;
                }

                push(u, a, min(excess[u], a.cap));
                if (excess[a.to] > 0 && !queued[a.to]) active.push(a.to), queued[a.to] = true;
            }
        }
    }

    long long reduced_cost(int u, const Arc& a) const{
        return potential[u] + a.cost - potential[a.to];
    }

    void push(int u, Arc& a, long long amount){
        a.cap -= amount, adj[a.to][a.rev].cap += amount;
        excess[u] -= amount, excess[a.to] += amount;
    }

    /// An excess node always has a residual arc: its excess arrived along arcs whose reverses are residual
    void relabel(int u, long long eps){
        long long best = LLONG_MIN;
        for (auto& a : adj[u]){
            if (a.cap > 0) best = max(best, potential[a.to] - a.cost);
        }
        potential[u] = best - eps;
        cur[u] = 0;
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

#define MAX 110

const long long INF = 1LL << 40, BIG = 1000;

int n, m, ex[MAX];

int main(){
    int i, u, v, c, f;
    long long base = 0, need = 0;

    read(n, m);
    MinCostCirculation g(n + 2);
    int S = n, T = n + 1;
    for (i = 0; i < m; i++){
        read(u, v, c, f);
        u--, v--;
        ex[v] += f, ex[u] -= f;
        if (f <= c){
            g.add_edge(u, v, c - f, 1);
            g.add_edge(u, v, INF, 2);
            g.add_edge(v, u, f, 1);
        }
        else{
            base += f - c;
            g.add_edge(u, v, INF, 2);
            g.add_edge(v, u, f - c, 0);
            g.add_edge(v, u, c, 1);
        }
    }
    g.add_edge(n - 1, 0, INF, 0);

    for (i = 0; i < n; i++){
        if (ex[i] > 0) g.add_edge(S, i, ex[i], -BIG), need += ex[i];
        if (ex[i] < 0) g.add_edge(i, T, -ex[i], 0);
    }
    g.add_edge(T, S, INF, 0);

    write(base + g.solve() + BIG * need);
    return 0;
}
