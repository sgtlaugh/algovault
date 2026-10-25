/***
 *
 * Steiner Tree (Dreyfus Wagner)
 * Minimum total weight of edges connecting every terminal of an undirected graph
 *
 * Complexity: O(3^k n + 2^k m log n) time, O(2^k n) memory for k terminals
 *
 * SteinerTree g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v, w): adds an undirected edge with w >= 0
 * g.solve(terminals): the minimum weight, 0 for fewer than two terminals, STEINER_INF if they are not all connected
 * g.solve_subsets(terminals): best, where best[mask] is solve of the terminals[i] with bit i set in mask
 *
 * Notes:
 *   - solve does not modify the graph, so one graph answers several terminal sets
 *   - Steiner forest, each group connected: forest[s] = min over t of forest[s ^ t] + best[t], t a union of groups in s
 *   - The total weight of all edges must stay below STEINER_INF (~2.3e18)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long STEINER_INF = numeric_limits<long long>::max() / 4;

struct SteinerTree{
    int n;
    vector<vector<pair<int, long long>>> adj;

    SteinerTree(int n) : n(n), adj(n) {}

    void add_edge(int u, int v, long long w){
        assert(w >= 0);
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    long long solve(const vector<int>& terminals) const{
        return solve_subsets(terminals).back();
    }

    vector<long long> solve_subsets(const vector<int>& terminals) const{
        int k = terminals.size();
        if (k <= 1) return vector<long long>(1 << k, 0);

        vector<vector<long long>> dp(1 << k, vector<long long>(n, STEINER_INF));
        for (int i = 0; i < k; i++) dp[1 << i][terminals[i]] = 0;

        for (int mask = 1; mask < (1 << k); mask++){
            vector<long long>& d = dp[mask];
            for (int sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask){
                if (sub < (sub ^ mask)) continue;  /// each split once
                for (int v = 0; v < n; v++) d[v] = min(d[v], dp[sub][v] + dp[sub ^ mask][v]);
            }

            priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> heap;
            for (int v = 0; v < n; v++){
                if (d[v] < STEINER_INF) heap.push({d[v], v});
            }
            while (!heap.empty()){
                long long dv = heap.top().first;
                int v = heap.top().second;
                heap.pop();
                if (dv != d[v]) continue;
                for (auto e : adj[v]){
                    int u = e.first;
                    long long w = e.second;
                    if (dv + w < d[u]) d[u] = dv + w, heap.push({d[u], u});
                }
            }
        }

        vector<long long> best(1 << k, 0);
        for (int mask = 1; mask < (1 << k); mask++) best[mask] = *min_element(dp[mask].begin(), dp[mask].end());
        return best;
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
    int t, n, m, k, i, s, sub, u, v, w;

    while (read(t)){
        while (t--){
            read(n, m, k);
            SteinerTree g(n);
            for (i = 0; i < m; i++){
                read(u, v, w);
                g.add_edge(u - 1, v - 1, w);
            }

            vector<int> terminals;
            for (i = 0; i < k; i++) terminals.push_back(i);
            for (i = 0; i < k; i++) terminals.push_back(n - k + i);
            vector<long long> best = g.solve_subsets(terminals);

            int full = (1 << (2 * k)) - 1, half = (1 << k) - 1;
            vector<long long> forest(full + 1, STEINER_INF);
            forest[0] = 0;
            for (s = 1; s <= full; s++){
                if (__builtin_popcount(s & half) != __builtin_popcount(s >> k)) continue;
                for (sub = s; sub > 0; sub = (sub - 1) & s){
                    if (__builtin_popcount(sub & half) != __builtin_popcount(sub >> k)) continue;
                    forest[s] = min(forest[s], forest[s ^ sub] + best[sub]);
                }
            }

            if (forest[full] >= STEINER_INF) write("No solution");
            else write(forest[full]);
        }
    }
    return 0;
}
