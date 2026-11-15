/***
 *
 * Monotonic Stack
 * Nearest smaller or greater element on each side, and the largest rectangle in a histogram
 *
 * Complexity: O(n)
 *
 * largest_rectangle(h): largest axis aligned rectangle under the histogram h, as long long
 * next_smaller(a)[i]: smallest j > i with a[j] < a[i], or n
 * next_smaller(a, compare)[i]: smallest j > i with compare(a[j], a[i]), or n
 * previous_smaller(a)[i]: largest j < i with a[j] < a[i], or -1
 * previous_smaller(a, compare)[i]: largest j < i with compare(a[j], a[i]), or -1
 *
 * Notes:
 *   - compare is strict like greater<T>() for the nearest greater element, or non-strict like less_equal<T>()
 *   - a[i] is the minimum of every range that contains i and lies inside (previous_smaller(a)[i], next_smaller(a)[i])
 *   - For the sum of subarray minimums, pass less_equal<T>() to one side so a tied minimum is counted once
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
vector<int> next_smaller(const vector<T>& a, Compare compare = Compare()){
    int n = a.size();
    vector<int> res(n), stack;

    for (int i = n - 1; i >= 0; i--){
        while (!stack.empty() && !compare(a[stack.back()], a[i])) stack.pop_back();
        res[i] = stack.empty() ? n : stack.back();
        stack.push_back(i);
    }

    return res;
}

template <typename T, typename Compare = less<T>>
vector<int> previous_smaller(const vector<T>& a, Compare compare = Compare()){
    int n = a.size();
    vector<int> res(n), stack;

    for (int i = 0; i < n; i++){
        while (!stack.empty() && !compare(a[stack.back()], a[i])) stack.pop_back();
        res[i] = stack.empty() ? -1 : stack.back();
        stack.push_back(i);
    }

    return res;
}

template <typename T>
long long largest_rectangle(const vector<T>& h){
    auto left = previous_smaller(h), right = next_smaller(h);
    long long best = 0;
    for (int i = 0; i < (int)h.size(); i++) best = max(best, (long long)h[i] * (right[i] - left[i] - 1));
    return best;
}

/***
 *
 * Sparse Table
 * O(1) range queries on a static array for operations that ignore repeats: min, max, gcd, and, or
 *
 * Complexity:
 *   - SparseTable: O(n log n) build time and memory, O(1) per query
 *   - LinearSparseTable: O(n) build time and memory, O(1) per query
 *
 * SparseTable<T> table(a): range minimum of a
 * SparseTable table(a, op): range op of a, for an op that ignores repeats like max or gcd
 * LinearSparseTable<T> table(a): range minimum of a
 * LinearSparseTable<T, greater<T>> table(a): range maximum of a
 * table.query(l, r): the answer over a[l..r], 0 <= l <= r < n
 *
 * Notes:
 *   - Sum, xor and product count repeats, use disjoint_sparse_table.cpp for them
 *   - At n = 1e6 LinearSparseTable builds 3x faster in 12 MB instead of 72 MB, and its queries are about 2x slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct MinOp{
    T operator()(const T& x, const T& y) const{
        return min(x, y);
    }
};

template <typename T, typename Op = MinOp<T>>
struct SparseTable{
    vector<vector<T>> dp;  /// dp[h][i] = op over a[i .. i + 2^h - 1]
    Op combine;

    SparseTable(const vector<T>& a, Op combine = Op()) : dp(1, a), combine(combine){
        int n = a.size();
        for (int h = 1; (1 << h) <= n; h++){
            dp.push_back(vector<T>(n - (1 << h) + 1));
            for (int i = 0; i + (1 << h) <= n; i++) dp[h][i] = combine(dp[h - 1][i], dp[h - 1][i + (1 << (h - 1))]);
        }
    }

    /// Two blocks of length 2^h cover [l, r], overlapping in the middle
    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)dp[0].size());
        int h = __lg(r - l + 1);
        return combine(dp[h][l], dp[h][r - (1 << h) + 1]);
    }
};

template <typename T, typename Compare = less<T>>
struct LinearSparseTable{
    static constexpr int B = 64;
    Compare compare;
    vector<T> val;
    vector<unsigned long long> mask;  /// bit k of mask[i]: val[i - k] beats every value after it up to i
    vector<vector<int>> table;        /// sparse table of block winners, by index

    LinearSparseTable(const vector<T>& a) : val(a), mask(a.size()){
        int n = val.size();
        unsigned long long stack = 0;  /// the monotonic stack of the last 64 positions, as a bitmask

        for (int i = 0; i < n; i++){
            stack <<= 1;
            /// the lowest set bit is the newest entry, pop while val[i] is at least as good
            while (stack && !compare(val[i - __builtin_ctzll(stack)], val[i])) stack &= stack - 1;
            stack |= 1;
            mask[i] = stack;
        }

        int blocks = n / B;  /// a query reads only blocks wholly before r's, so a partial last block is never read
        table.assign(1, vector<int>(blocks));
        for (int b = 0; b < blocks; b++) table[0][b] = best_ending_at(b * B + B - 1);
        for (int k = 1; (1 << k) <= blocks; k++){
            table.push_back(vector<int>(blocks - (1 << k) + 1));
            for (int b = 0; b + (1 << k) <= blocks; b++) table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
        }
    }

    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)val.size());
        if (r - l + 1 <= B) return val[best_ending_at(r, r - l + 1)];

        int res = better(best_ending_at(l + B - 1), best_ending_at(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = __lg(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }

        return val[res];
    }

    /// Index of the best value in [r - len + 1, r], len <= 64: the oldest stack entry inside the window
    int best_ending_at(int r, int len = B) const{
        unsigned long long m = len == B ? mask[r] : mask[r] & ((1ULL << len) - 1);
        return r - __lg(m);
    }

    int better(int i, int j) const{
        return compare(val[j], val[i]) ? j : i;
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

const int ULTRA = 150000;

int main(){
    int n, i, x;

    while (read(n)){
        vector<int> h(n);
        read(h);

        vector<int> left = previous_smaller(h, greater<int>()), right = next_smaller(h, greater<int>());
        SparseTable<int> table(h);

        vector<int> res;
        for (i = 0; i < n; i++){
            x = 0;
            if (left[i] != -1) x = max(x, table.query(left[i], i));
            if (right[i] != n) x = max(x, table.query(i, right[i]));
            if ((h[i] - x) >= ULTRA) res.push_back(i + 1);
        }

        write(res);
    }

    return 0;
}
