/***
 *
 * Floor Sum
 * Sums, minima and range counts over the arithmetic progression (a * i + b) mod m, by Euclid-like reduction
 *
 * Complexity: O(log m) per call
 *
 * count_mod_in_range(n, m, a, l, r): number of 0 <= x < n with l <= (a * x) mod m <= r, any n >= 0
 * first_mod_in_range(m, a, l, r): smallest x >= 0 with l <= (a * x) mod m <= r, or -1 if there is none
 * floor_sum(n, m, a, b): sum of floor((a * i + b) / m) for 0 <= i < n, n >= 0, m >= 1, any signed a, b
 *     Exact for any long long a, b when n <= 2^32, and for any n >= 0 when 0 <= a < m
 *     Returns __int128, cast to long long when the sum is known to fit
 * min_mod_linear(n, m, a, b): min of (a * i + b) mod m over 0 <= i < n, n >= 1, any signed a, b
 *
 * Notes:
 *   - Arguments are long long, l and r may lie outside [0, m - 1], they are clamped to it
 *   - Count of 0 <= x < n, y >= 0 with a x + b y <= c is floor_sum(n, b, -a, c + b) if a >= 0 and a (n - 1) <= c + b
 *   - Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

typedef unsigned __int128 u128;

/// Called with 0 <= a, b < m, so a * n + b < m * (n + 1) and the sum stays below n^2: 64 bits suffice for n, m < 2^32
template<typename U>
U floor_sum_unsigned(U n, U m, U a, U b){
    U res = 0;

    while (true){
        if (a >= m){
            res += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m){
            res += n * (b / m);
            b %= m;
        }

        U y_max = a * n + b;
        if (y_max < m) return res;
        n = y_max / m, b = y_max % m;
        swap(m, a);
    }
}

long long floor_sum_reduce_mod(long long x, long long m){
    x %= m;
    return x < 0 ? x + m : x;
}

__int128 floor_sum(long long n, long long m, long long a, long long b){
    assert(n >= 0 && m >= 1);

    long long ra = floor_sum_reduce_mod(a, m), rb = floor_sum_reduce_mod(b, m);
    long long qa = a / m - (a % m < 0), qb = b / m - (b % m < 0);
    __int128 res = (__int128)qa * ((__int128)n * (n - 1) / 2) + (__int128)qb * n;

    /// The 64-bit path is 46% faster for n, m <= 1e9 and cannot overflow below 2^32, see floor_sum_unsigned
    if (n < (1LL << 32) && m < (1LL << 32)) return res + floor_sum_unsigned<unsigned long long>(n, m, ra, rb);
    return res + (__int128)floor_sum_unsigned<u128>(n, m, ra, rb);
}

long long count_mod_in_range(long long n, long long m, long long a, long long l, long long r){
    assert(n >= 0 && m >= 1);
    l = max(l, 0LL), r = min(r, m - 1);
    if (l > r || !n) return 0;

    a = floor_sum_reduce_mod(a, m);
    return floor_sum(n, m, a, m - l) - floor_sum(n, m, a, m - 1 - r);
}

/// Requires 0 <= a < m and 0 <= l <= r < m
long long first_mod_in_range_reduced(long long m, long long a, long long l, long long r){
    if (!a) return l ? -1 : 0;

    __int128 c = ((__int128)l + a - 1) / a;
    if (a * c <= r) return c;

    /// No multiple of a lies in [l, r]; with m = k * a + b, a * x - m * y in [l, r] gives -r <= (b * y) mod a <= -l
    long long b = m % a;
    long long y = first_mod_in_range_reduced(a, b, a - r % a, a - l % a);
    if (y == -1) return -1;
    return ((__int128)b * y + l + a - 1) / a + (__int128)(m / a) * y;
}

long long first_mod_in_range(long long m, long long a, long long l, long long r){
    assert(m >= 1);
    l = max(l, 0LL), r = min(r, m - 1);
    if (l > r) return -1;

    a = floor_sum_reduce_mod(a, m);
    return first_mod_in_range_reduced(m, a, l, r);
}

/// Each round keeps only the candidates that can be minimal, a progression mod a when 2a <= m, else mod m - a
long long min_mod_linear(long long n, long long m, long long a, long long b){
    assert(n >= 1 && m >= 1);
    a = floor_sum_reduce_mod(a, m);
    b = floor_sum_reduce_mod(b, m);

    long long res = m;
    while (a){
        if (2 * (__int128)a <= m){
            res = min(res, b);
            long long wraps = ((__int128)a * (n - 1) + b) / m;
            if (!wraps) return res;

            n = wraps;
            b = ((b - m) % a + a) % a;
            long long step = (a - m % a) % a;
            m = a, a = step;
        }
        else{
            long long c = m - a;
            res = min(res, (long long)(((__int128)a * (n - 1) + b) % m));
            if (b >= (__int128)c * n) return res;

            n = ((__int128)c * n - 1 - b) / m + 1;
            a = m % c, b %= c, m = c;
        }
    }

    return min(res, b);
}

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
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
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

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
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

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
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

const int coef[] = {11, 9, 81};

long long count_solutions(long long total){
    if (total < 0) return 0;

    long long a = (5 * total) % 9;
    if ((11 * a) > total) return 0;

    long long k = (total - 11 * a) / 9, steps = k / 11;
    return (long long)floor_sum(steps + 1, 9, 11, k - 11 * steps) + steps + 1;
}

long long solve(long long p, long long q, long long L){
    if (729 % q) return 0;

    long long total = 729 / q * p, cur, res = 0;
    for (int mask = 0; mask < 8; mask++){
        cur = total;
        for (int i = 0; i < 3; i++){
            if ((mask >> i) & 1) cur -= coef[i] * (L + 1);
        }

        if (__builtin_popcount(mask) & 1) res -= count_solutions(cur);
        else res += count_solutions(cur);
    }

    return res;
}

int main(){
    int t, T = 0;
    long long p, q, L;
    char buf[64];

    while (read(t)){
        while (t--){
            read(p, q, L);
            snprintf(buf, sizeof(buf), "Case %d: %lld", ++T, solve(p, q, L));
            write(string(buf));
        }
    }

    return 0;
}
