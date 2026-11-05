/***
 *
 * Plug DP
 * Counts or optimizes cycles and paths on a grid: the caller writes one cell's transition, the engine runs the rest
 *
 * Complexity: O(S) expected per cell for S states alive at that cell, plus the transition's own cost
 *
 * PlugDP<T> dp(n, m, bits, start): n x m grid, bits per plug, the all-zero state holds start, equal states add up
 * PlugDP<T, Merge> dp(n, m, bits, start, merge): equal states combine their values with merge(a, b), such as min
 * dp.add(s, v): only inside the transition, offers state s with value v after the current cell
 * dp.run(transition): calls transition(row, col, s, v) for each state s of value v before each cell, row-major, once
 * dp.get(s, pos): plug pos of state s
 * dp.join(s): s with the two bracket plugs entering the current cell cleared and their paths joined into one
 * dp.left(s): plug of state s entering the current cell from the left
 * dp.match(s, pos): position of the bracket paired with plug pos, plugs 1 opening and 2 closing, balanced
 * dp.normalize(s): s with its nonzero labels renumbered 1, 2, 3... in order of first appearance from position 0
 * dp.place(s, down, right): s with the plugs entering the current cell replaced by down and right, both below 2^bits
 * dp.result(missing): EDGES only, value of the all-zero state after the last cell, missing (T() by default) if none
 * dp.set(s, pos, x): state s with plug pos replaced by x, 0 <= x < 2^bits
 * dp.up(s): plug of state s entering the current cell from above
 * dp.profile: EDGES by default, set dp.profile = dp.CELLS before run when each plug holds a cell's label
 * dp.states, dp.values: the states alive before the current cell, or after the last one, and their values
 *
 * Notes:
 *   - Before cell (row, col) plug col is its edge from the left and col + 1 from above, after it down and right
 *   - left, up, place and join read the current cell, so like add they work only inside the transition
 *   - join asserts both plugs are 1 opening or 2 closing, never left 1 with up 2: that closes a cycle, place(s, 0, 0)
 *   - join keeps brackets balanced: on left 1, up 1 the up plug's partner becomes 1, on 2, 2 the left plug's becomes 2
 *   - EDGES: add drops states with a nonzero plug leaving the grid, so transitions never test the border
 *   - CELLS: a cell's label rides on its down plug, place(s, label, 0), so before a cell up(s) is the cell above
 *   - CELLS: the cell to the left is get(s, col - 1), plug col is 0, and nothing is dropped at the border
 *   - CELLS: after run plug j of each state is the label of cell (n - 1, j), so read states and values, not result
 *   - Connectivity like spanning trees: cell labels keep far fewer states, 1430 against 30804 per row at m = 8
 *   - Transpose so m is the narrower side: UVa 10572's empty 3 x 8 grid visits 295 states at m = 3, 22319 at m = 8
 *   - A transition reading only col, s and the cell's kind can be memoized per column by the caller: 10x on UVa 10572
 *   - Small states hash too: a direct-indexed table measured only 1.1-1.2x faster at up to 20 state bits
 *   - (m + 1) bits <= 64, bits <= 8: 2 bits fit brackets, labels need 2^bits >= components + 2 to leave a fresh one
 *   - A state holds plugs 0..m only: the row shift moves every bit up by bits, so nothing can be kept above plug m
 *   - main() writes Formula 1 (one cycle, bracket encoding) and HDU 1693 (any cycles, 1 bit) as transitions
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Merge = plus<T>>
struct PlugDP{
    enum Profile{EDGES, CELLS};

    int n, m, bits, row = 0, col = 0, lg = 10;
    uint64_t off_grid = 0;
    uint64_t low = 0;  /// bit 0 of plugs 0..m, so x * low holds x in every plug
    Merge merge;
    Profile profile = EDGES;
    vector<uint64_t> states, next_states;
    vector<T> values, next_values;
    vector<int> head, link;

    PlugDP(int n, int m, int bits, T start, Merge merge = Merge())
        : n(n), m(m), bits(bits), merge(merge), states{0}, values{start}, head(1 << lg, -1){
        assert(n >= 1 && m >= 1 && bits >= 1 && bits <= 8 && (m + 1) * bits <= 64);
        for (int p = 0; p <= m; p++) low |= uint64_t(1) << (p * bits);
    }

    void add(uint64_t s, T v){
        if (s & off_grid) return;

        size_t h = slot(s);
        for (int i = head[h]; i >= 0; i = link[i]){
            if (next_states[i] == s){
                next_values[i] = merge(next_values[i], v);
                return;
            }
        }

        link.push_back(head[h]);
        head[h] = next_states.size();
        next_states.push_back(s);
        next_values.push_back(v);
    }

    template <typename Transition>
    void run(Transition transition){
        uint64_t ones = profile == EDGES ? (uint64_t(1) << bits) - 1 : 0;
        for (row = 0; row < n; row++){
            for (col = 0; col < m; col++){
                /// sized per cell, not per insert: a resize check inside add stops GCC inlining it, 1.6x slower
                while (head.size() < 8 * states.size()) head.assign(size_t(1) << ++lg, -1);
                off_grid = (col == m - 1 ? ones << (m * bits) : 0) | (row == n - 1 ? ones << (col * bits) : 0);
                for (size_t i = 0; i < states.size(); i++) transition(row, col, states[i], values[i]);

                for (uint64_t s : next_states) head[slot(s)] = -1;
                link.clear();
                swap(states, next_states), swap(values, next_values);
                next_states.clear(), next_values.clear();
            }
            /// not after the last row, so a CELLS state ends with cell (n - 1, j) at plug j
            if (row + 1 < n){
                for (uint64_t& s : states) s <<= bits;
            }
        }
    }

    int get(uint64_t s, int pos) const{
        return s >> (pos * bits) & ((1 << bits) - 1);
    }

    uint64_t join(uint64_t s) const{
        int l = left(s), u = up(s);
        assert(l >= 1 && l <= 2 && u >= 1 && u <= 2 && !(l == 1 && u == 2));

        uint64_t t = place(s, 0, 0);
        if (l == 1 && u == 1) return set(t, match(s, col + 1), 1);
        if (l == 2 && u == 2) return set(t, match(s, col), 2);
        return t;
    }

    int left(uint64_t s) const{
        return get(s, col);
    }

    int match(uint64_t s, int pos) const{
        /// a shift per plug, not get's multiply: 1.4x faster on random brackets; s == 0 leaves no partner to find
        int mask = (1 << bits) - 1, depth = 0;
        if (get(s, pos) == 1){
            for (s >>= pos * bits; s; s >>= bits, pos++){
                int x = s & mask;
                depth += (x == 1) - (x == 2);
                if (depth == 0) return pos;
            }
        }
        else{
            for (int shift = pos * bits; shift >= 0; shift -= bits, pos--){
                int x = s >> shift & mask;
                depth += (x == 1) - (x == 2);
                if (depth == 0) return pos;
            }
        }

        assert(false);
        return -1;
    }

    uint64_t normalize(uint64_t s) const{
        /// a step per label, not per plug, 1.6x faster; a fresh word, so a new id never meets a label not yet renamed
        uint64_t left = nonzero(s), res = 0;
        for (uint64_t id = 1; left; id++){
            uint64_t same = equal(s, s >> __builtin_ctzll(left) & ((1 << bits) - 1));
            res |= same * id;
            left &= ~same;
        }

        return res;
    }

    uint64_t place(uint64_t s, int down, int right) const{
        /// one mask over both plugs instead of two set calls: 4% faster on the HDU 1693 transition, measured
        int shift = col * bits;
        uint64_t both = (uint64_t(1) << 2 * bits) - 1;
        return (s & ~(both << shift)) | (uint64_t(down) | uint64_t(right) << bits) << shift;
    }

    T result(T missing = T()) const{
        assert(profile == EDGES);
        return states.empty() ? missing : values[0];
    }

    uint64_t set(uint64_t s, int pos, int x) const{
        int shift = pos * bits;
        return (s & ~(uint64_t((1 << bits) - 1) << shift)) | uint64_t(x) << shift;
    }

    int up(uint64_t s) const{
        return get(s, col + 1);
    }

    size_t slot(uint64_t s) const{
        return s * 0x9E3779B97F4A7C15ULL >> (64 - lg);
    }

    uint64_t equal(uint64_t s, uint64_t x) const{
        return low & ~nonzero(s ^ x * low);
    }

    uint64_t nonzero(uint64_t s) const{
        uint64_t x = s;
        for (int b = 1; b < bits; b++) x |= s >> b;
        return x & low;
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

int n, m, last;
vector<string> grid;

void transpose(){
    vector<string> res(m, string(n, '.'));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) res[j][i] = grid[i][j];
    }

    grid = res;
    swap(n, m);
}

int main(){
    int i;

    read(n, m);
    grid.resize(n);
    read(grid);
    if (m > n) transpose();

    last = -1;
    for (i = 0; i < n * m; i++){
        if (grid[i / m][i % m] == '.') last = i;
    }

    PlugDP<long long> dp(n, m, 2, 1);
    dp.run([&](int r, int c, uint64_t s, long long v){
        int l = dp.left(s), u = dp.up(s);
        if (grid[r][c] == '*'){
            if (!l && !u) dp.add(s, v);
        }
        else if (!l && !u) dp.add(dp.place(s, 1, 2), v);
        else if (!l || !u){
            dp.add(dp.place(s, l | u, 0), v);
            dp.add(dp.place(s, 0, l | u), v);
        }
        else if (l == 1 && u == 2){
            if ((r * m + c) == last) dp.add(dp.place(s, 0, 0), v);
        }
        else dp.add(dp.join(s), v);
    });

    write(dp.result());
    return 0;
}
