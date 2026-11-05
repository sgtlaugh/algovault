/***
 *
 * 2SAT with Kosaraju's algorithm
 * Finds a true or false value for every variable that satisfies a system of constraints on pairs of variables
 *
 * Complexity: O(n + m), n = number of variables, m = number of clauses
 *
 * TwoSat g(n): no constraints on variables 1 to n, literal x means x is true and -x means x is false
 * g.add_and(a, b): requires both literals a and b
 * g.add_implication(a, b): requires b whenever a holds
 * g.add_or(a, b): requires at least one of a and b
 * g.add_var(): adds variable n + 1 and returns its index
 * g.add_xor(a, b): requires exactly one of a and b
 * g.at_most_one(literals): allows at most one of the literals to hold, a literal listed twice counts twice
 * g.force_false(x): requires literal x to be false
 * g.force_true(x): requires literal x to be true
 * g.get_assignment(): the true variables, ascending, valid after is_satisfiable() returns true
 * g.is_satisfiable(): whether the constraints can all hold, finding an assignment when they can
 * g.value(x): whether literal x holds in that assignment, valid after is_satisfiable() returns true
 *
 * Notes:
 *   - at_most_one adds k - 2 variables via add_var() and 3k - 5 clauses for k >= 2, not k^2 / 2 pairwise clauses
 *   - Variables from add_var() are numbered after the existing ones, so earlier indices never change
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TwoSat{
    int n, t;
    vector<pair<int, int>> edges;  /// signed literals, mapped to nodes 1..2n only in build() since add_var() changes n
    vector<int> adj_start, adj_to, rev_start, rev_to;
    vector<int> visited, comp, order, dfs_t;

    TwoSat(int n): n(n) {}

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    inline void add_implication(int a, int b){
        assert(a != 0 && abs(a) <= n && b != 0 && abs(b) <= n);
        edges.push_back({a, b}), edges.push_back({-b, -a});
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    int add_var(){
        return ++n;
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    /// prev is a literal meaning "one of lits[0..i-1] is true", it forbids lits[i] and carries over to the next prefix
    void at_most_one(const vector<int>& lits){
        int k = lits.size();

        if (k <= 1) return;

        int prev = lits[0];
        for (int i = 1; i + 1 < k; i++){
            int cur = add_var();
            add_implication(prev, -lits[i]);
            add_implication(prev, cur);
            add_implication(lits[i], cur);
            prev = cur;
        }

        add_implication(prev, -lits[k - 1]);
    }

    inline void force_false(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({x, -x});  /// x -> !x is its own contrapositive, add it once
    }

    inline void force_true(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({-x, x});  /// !x -> x is its own contrapositive, add it once
    }

    vector<int> get_assignment(){
        vector<int> res;
        for (int x = 1; x <= n; x++){
            if (value(x)) res.push_back(x);
        }

        return res;
    }

    bool is_satisfiable(){
        build();
        for (int i = 1; i <= n; i++){
            if (comp[i] == comp[i + n]) return false;
        }

        return true;
    }

    /// If x implies !x then comp[!x] is discovered first, so x is false
    inline bool value(int x){
        if (x < 0) return !value(-x);
        return comp[x] < comp[x + n];
    }

    /// Components are numbered in discovery order, which is a reverse topological order of the implication graph
    void build(){
        int m = 2 * n + 2, c = 0;

        flatten(adj_start, adj_to, false);
        flatten(rev_start, rev_to, true);
        dfs_t.assign(m, 0), order.assign(m, 0), comp.assign(m, 0);

        vector<pair<int, int>> topsort_stack;
        visited.assign(m, 0);
        t = 0;
        for (int i = 2 * n; i >= 1; i--){
            if (!visited[i]) topsort(i, topsort_stack);
            order[dfs_t[i]] = i;
        }

        vector<int> dfs_stack;
        visited.assign(m, 0);
        for (int i = 2 * n; i >= 1; i--){
            int x = order[i];
            if (!visited[x]) dfs(x, c++, dfs_stack);
        }
    }

    void flatten(vector<int>& start, vector<int>& to, bool reverse){
        start.assign(2 * n + 2, 0), to.resize(edges.size());
        for (auto [a, b]: edges) start[node(reverse ? b : a)]++;
        for (int i = 1; i <= 2 * n + 1; i++) start[i] += start[i - 1];
        for (auto [a, b]: edges){
            int u = node(a), v = node(b);
            if (reverse) swap(u, v);
            to[--start[u]] = v;
        }
    }

    /// Literal x is node x, literal -x is node n + x
    inline int node(int x){
        return x > 0 ? x : n - x;
    }

    void topsort(int s, vector<pair<int, int>>& stack){
        visited[s] = true, stack.push_back({s, rev_start[s]});
        while (!stack.empty()){
            auto& [u, e] = stack.back();
            if (e < rev_start[u + 1]){
                int v = rev_to[e++];
                if (!visited[v]) visited[v] = true, stack.push_back({v, rev_start[v]});
            }
            else{
                dfs_t[u] = ++t;
                stack.pop_back();
            }
        }
    }

    /// Labels everything reachable from s, the visiting order does not matter here
    void dfs(int s, int c, vector<int>& stack){
        comp[s] = c, visited[s] = true, stack.push_back(s);
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            for (int e = adj_start[u]; e < adj_start[u + 1]; e++){
                int v = adj_to[e];
                if (!visited[v]) comp[v] = c, visited[v] = true, stack.push_back(v);
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

#define MAX 110

const int dx[] = {0, 0, -1, 1};
const int dy[] = {-1, 1, 0, 0};

int r, c, id[MAX][MAX];
vector<string> grid;

int scan(int i, int j, int d){
    for (i += dx[d], j += dy[d]; i >= 0 && i < r && j >= 0 && j < c; i += dx[d], j += dy[d]){
        if (grid[i][j] == '#') return 0;
        if (grid[i][j] == 'T') return id[i][j];
    }

    return 0;
}

int main(){
    int i, j, k, x, y;

    read(r, c);
    grid.resize(r);
    read(grid);

    k = 0;
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            if (grid[i][j] == 'T') id[i][j] = ++k;
        }
    }

    TwoSat sat(2 * k);
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            if (grid[i][j] == 'T'){
                x = id[i][j];
                if (scan(i, j, 0)) sat.force_true(x);
                if (scan(i, j, 1)) sat.force_false(x);
                if (scan(i, j, 2)) sat.force_false(x + k);
                if (scan(i, j, 3)) sat.force_true(x + k);
            }
            if (grid[i][j] != 'n') continue;

            vector<int> lits;
            x = scan(i, j, 0), y = scan(i, j, 1);
            if (x && !y) lits.push_back(x);
            if (y && !x) lits.push_back(-y);
            x = scan(i, j, 2), y = scan(i, j, 3);
            if (x && !y) lits.push_back(-(x + k));
            if (y && !x) lits.push_back(y + k);

            if (lits.size() == 1) sat.force_true(lits[0]);
            else sat.add_or(lits[0], lits[1]);
        }
    }

    sat.is_satisfiable();
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            if (grid[i][j] != 'T') continue;
            x = id[i][j];
            if (sat.value(x)) grid[i][j] = sat.value(x + k) ? '3' : '2';
            else grid[i][j] = sat.value(x + k) ? '4' : '1';
        }
        write(grid[i]);
    }

    return 0;
}
