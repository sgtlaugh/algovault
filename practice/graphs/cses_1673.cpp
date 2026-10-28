/***
 *
 * Bellman Ford
 * Single source shortest paths on a directed graph with negative edges, and negative cycle retrieval
 *
 * Complexity: O(n m) per call
 *
 * BellmanFord g(n): nodes are 0 to n - 1
 * g.add_edge(u, v, w): directed edge u -> v with long long weight w, negative allowed
 * g.negative_cycle(): one negative cycle anywhere in the graph, cycle[i] -> cycle[i + 1] an edge, empty if none
 * g.shortest_paths(src)[v]: distance from src, INF if unreachable, NEG_INF if a negative cycle on the way reaches v
 *
 * Notes:
 *   - Requires n max|w| <= 9e18
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct BellmanFord{
    struct Edge{
        int u, v;
        long long w;
    };

    static constexpr long long INF = numeric_limits<long long>::max();
    static constexpr long long NEG_INF = numeric_limits<long long>::min();

    int n;
    vector<Edge> edges;

    BellmanFord(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    vector<int> negative_cycle() const{
        long long low = -(n - 1) * max_abs_weight();
        vector<long long> dist(n, 0);
        vector<int> parent(n, -1);

        for (int pass = 0; pass < n; pass++){
            int relaxed = -1;
            for (auto [u, v, w] : edges){
                if (dist[u] + w >= dist[v]) continue;
                dist[v] = dist[u] + w, parent[v] = u, relaxed = v;
                /// Below every simple path weight, so the parent chain of v cannot end at a root: it holds a cycle
                if (dist[v] < low) return cycle_through(parent, v);
            }
            if (relaxed == -1) return {};
            /// An acyclic parent chain has at most n - 1 edges, so a relaxation in pass n means the chain holds a cycle
            if (pass == n - 1) return cycle_through(parent, relaxed);
        }

        return {};
    }

    vector<long long> shortest_paths(int src) const{
        long long low = -(n - 1) * max_abs_weight();
        vector<long long> dist(n, INF);
        dist[src] = 0;

        /// Unaffected vertices are exact after n - 1 passes, an edge relaxing after that is in a negative cycle's reach
        for (int pass = 0; pass < n; pass++){
            bool changed = false;
            for (auto [u, v, w] : edges){
                if (dist[u] == INF || dist[u] == NEG_INF || dist[u] + w >= dist[v]) continue;
                dist[v] = (pass == n - 1 || dist[u] + w < low) ? NEG_INF : dist[u] + w;
                changed = true;
            }
            if (!changed) break;
        }

        for (bool changed = true; changed;){
            changed = false;
            for (auto [u, v, w] : edges){
                if (dist[u] == NEG_INF && dist[v] != NEG_INF) dist[v] = NEG_INF, changed = true;
            }
        }

        return dist;
    }

    long long max_abs_weight() const{
        long long res = 0;
        for (auto [u, v, w] : edges) res = max(res, abs(w));
        return res;
    }

    vector<int> cycle_through(const vector<int>& parent, int v) const{
        for (int i = 0; i < n; i++) v = parent[v];

        vector<int> cycle = {v};
        for (int u = parent[v]; u != v; u = parent[u]) cycle.push_back(u);
        reverse(cycle.begin(), cycle.end());
        return cycle;
    }
};

bool is_rotation(vector<int> cycle, const vector<int>& expected){
    for (size_t i = 0; i < cycle.size(); i++){
        if (cycle == expected) return true;
        rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
    }

    return false;
}

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
    int n, m, i, a, b, x;

    read(n, m);
    BellmanFord g(n);
    for (i = 0; i < m; i++){
        read(a, b, x);
        g.add_edge(a - 1, b - 1, -x);
    }

    long long res = g.shortest_paths(0)[n - 1];
    write(res == BellmanFord::NEG_INF ? -1 : -res);
    return 0;
}
