/***
 *
 * Minimum Path Cover in a DAG
 * Vertex disjoint path cover, path cover with shared vertices, and maximum antichain (Dilworth)
 *
 * Complexity:
 *   - min_disjoint_path_cover: O(m sqrt(n))
 *   - min_path_cover and max_antichain: O(n^3 / 64 + n^2.5), O(n^2) memory
 *
 * DAGPathCover g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v): directed edge, the graph must stay acyclic
 * g.max_antichain(): a largest set of pairwise unreachable vertices, its size equals min_path_cover()
 * g.min_disjoint_path_cover(): fewest paths covering every vertex exactly once
 * g.min_path_cover(): fewest paths covering every vertex at least once, paths may share vertices
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY hopcroft_karp from code_library/graphs/hopcroft_karp.cpp
struct HopcroftKarp{
    int n_left, n_right;
    vector<vector<int>> adj;
    vector<int> match_left, match_right, dist, it;

    HopcroftKarp(int n_left, int n_right) : n_left(n_left), n_right(n_right), adj(n_left), match_left(n_left, -1), match_right(n_right, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    int max_matching(){
        while (true){
            int limit = bfs();
            if (limit == INT_MAX) break;
            it.assign(n_left, 0);
            for (int u = 0; u < n_left; u++){
                if (match_left[u] == -1) augment(u, limit);
            }
        }

        return n_left - count(match_left.begin(), match_left.end(), -1);
    }

    /// Layers left vertices by distance from free left ones, stopping at the first layer reaching a free right vertex
    int bfs(){
        dist.assign(n_left, -1);
        vector<int> queue;
        for (int u = 0; u < n_left; u++){
            if (match_left[u] == -1) dist[u] = 0, queue.push_back(u);
        }

        int limit = INT_MAX;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int v : adj[u]){
                int w = match_right[v];
                if (w == -1) limit = min(limit, dist[u]);
                else if (dist[w] == -1 && dist[u] < limit){
                    dist[w] = dist[u] + 1;
                    queue.push_back(w);
                }
            }
        }

        return limit;
    }

    void augment(int root, int limit){
        vector<int> stack = {root};
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == (int)adj[u].size()){
                dist[u] = -1;
                stack.pop_back();
                continue;
            }

            int v = adj[u][it[u]++], w = match_right[v];
            if (w == -1 && dist[u] == limit){
                for (int x : stack){
                    int y = adj[x][it[x] - 1];
                    match_left[x] = y, match_right[y] = x;
                }

                return;
            }
            if (w != -1 && dist[w] == dist[u] + 1) stack.push_back(w);
        }
    }
};
/// END COPY hopcroft_karp

struct DAGPathCover{
    int n;
    vector<vector<int>> adj;

    DAGPathCover(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    vector<int> max_antichain() const{
        HopcroftKarp hk = matching(n, closure());

        /// Konig: alternating reachability from free left vertices gives the minimum vertex cover
        vector<char> seen_left(n, 0), seen_right(n, 0);
        vector<int> queue;
        for (int u = 0; u < n; u++){
            if (hk.match_left[u] == -1) seen_left[u] = 1, queue.push_back(u);
        }
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : hk.adj[queue[i]]){
                if (seen_right[v]) continue;
                seen_right[v] = 1;
                int w = hk.match_right[v];
                if (w != -1 && !seen_left[w]) seen_left[w] = 1, queue.push_back(w);
            }
        }

        /// The cover is unvisited left plus visited right copies, the antichain is every vertex with neither copy in it
        vector<int> res;
        for (int v = 0; v < n; v++){
            if (seen_left[v] && !seen_right[v]) res.push_back(v);
        }

        return res;
    }

    int min_disjoint_path_cover() const{
        return n - matching(n, adj).max_matching();
    }

    int min_path_cover() const{
        return n - matching(n, closure()).max_matching();
    }

    /// Split every vertex into an out copy (left) and an in copy (right), each matched edge joins two paths
    static HopcroftKarp matching(int n, vector<vector<int>> edges){
        HopcroftKarp hk(n, n);
        hk.adj = move(edges);
        hk.max_matching();
        return hk;
    }

    vector<vector<int>> closure() const{
        int words = (n + 63) / 64;
        vector<vector<unsigned long long>> reach(n, vector<unsigned long long>(words, 0));

        /// Reverse topological order, so every successor's reach set is complete before it is merged
        vector<int> indegree(n, 0), order;
        for (int u = 0; u < n; u++){
            for (int v : adj[u]) indegree[v]++;
        }
        for (int u = 0; u < n; u++){
            if (!indegree[u]) order.push_back(u);
        }
        for (int i = 0; i < (int)order.size(); i++){
            for (int v : adj[order[i]]){
                if (!--indegree[v]) order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i >= 0; i--){
            int u = order[i];
            for (int v : adj[u]){
                reach[u][v >> 6] |= 1ULL << (v & 63);
                for (int w = 0; w < words; w++) reach[u][w] |= reach[v][w];
            }
        }

        vector<vector<int>> res(n);
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++){
                if (reach[u][v >> 6] >> (v & 63) & 1) res[u].push_back(v);
            }
        }

        return res;
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

int main(){
    int n, i, j, x;

    read(n);
    vector<int> ar;
    for (i = 0; i < n; i++){
        read(x);
        if (find(ar.begin(), ar.end(), x) == ar.end()) ar.push_back(x);
    }

    n = ar.size();
    DAGPathCover g(n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i != j && ar[j] % ar[i] == 0) g.add_edge(i, j);
        }
    }

    vector<int> res;
    for (int v: g.max_antichain()) res.push_back(ar[v]);
    write((int)res.size());
    write(res);
    return 0;
}
