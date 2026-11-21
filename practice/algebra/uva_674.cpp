/***
 *
 * Sparse Power Series
 * Inverse, log, exp, power and square root of a power series with few nonzero terms, modulo any prime
 *
 * Complexity: O(n k) for the first n coefficients of a series with k nonzero terms, plus O(log mod) per call
 *
 * SparseSeries f(terms, mod = 998244353): the sum of c x^d over the pairs {d, c} of terms, repeated degrees add up
 * f.exp(n): e^f mod x^n, requires a zero constant term
 * f.inverse(n): 1 / f mod x^n, requires a nonzero constant term
 * f.log(n): ln f mod x^n, requires the constant term 1
 * f.pow(k, n): f^k mod x^n for any k >= 0, f may start with zeros, and 0^0 = 1
 * f.sqrt(n): a square root of f mod x^n, lowest nonzero coefficient the smaller modular root, nullopt when none exists
 * f.terms: the nonzero terms {d, c} by increasing degree, c in [1, mod)
 *
 * Notes:
 *   - mod is a prime below 2^31 with n <= mod, so 1..n - 1 have inverses; sqrt also needs mod > 2
 *   - Degrees are >= 0 and coefficients any long long; results have exactly n coefficients in [0, mod)
 *   - Every call reads only the terms of degree below n
 *   - Modulo an NTT prime, polynomial.cpp does the same in O(n log n), faster once k passes about 80 at n = 10^6
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SparseSeries{
    vector<pair<int, int>> terms;
    int mod;

    SparseSeries(const vector<pair<int, long long>>& raw, int mod = 998244353) : mod(mod){
        assert(mod >= 2);
        vector<pair<int, long long>> sorted = raw;
        sort(sorted.begin(), sorted.end());

        for (auto& term : sorted){
            int d = term.first;
            long long c = term.second;
            assert(d >= 0);
            int r = (c % mod + mod) % mod;
            if (!terms.empty() && terms.back().first == d) terms.back().second = (terms.back().second + (long long)r) % mod;
            else terms.push_back({d, r});
        }
        terms.erase(remove_if(terms.begin(), terms.end(), [](const pair<int, int>& t){ return t.second == 0; }), terms.end());
    }

    vector<int> exp(int n) const{
        assert(n >= 0 && (terms.empty() || terms[0].first > 0));
        vector<int> res(n), inv = inverses(n);
        if (n == 0) return res;

        /// From g' = f' g: i g[i] = sum of d c g[i - d]
        res[0] = 1;
        for (int i = 1; i < n; i++){
            long long sum = 0;
            for (auto& term : terms){
                int d = term.first, c = term.second;
                if (d > i) break;
                sum = (sum + (long long)d * c % mod * res[i - d]) % mod;
            }
            res[i] = sum * inv[i] % mod;
        }

        return res;
    }

    vector<int> inverse(int n) const{
        assert(n >= 0 && !terms.empty() && terms[0].first == 0);
        vector<int> res(n);
        if (n == 0) return res;

        long long c0_inv = power(terms[0].second, mod - 2);
        res[0] = c0_inv;
        for (int i = 1; i < n; i++){
            long long sum = 0;
            for (int j = 1; j < (int)terms.size() && terms[j].first <= i; j++){
                sum = (sum + (long long)terms[j].second * res[i - terms[j].first]) % mod;
            }
            res[i] = (mod - sum) * c0_inv % mod;
        }

        return res;
    }

    vector<int> log(int n) const{
        assert(n >= 0 && !terms.empty() && terms[0] == make_pair(0, 1));
        vector<int> res(n), inv = inverses(n);

        /// q = (ln f)' solves f q = f', so q[m] = (m + 1) f[m + 1] - sum of c q[m - d] over d > 0
        vector<long long> q(max(n - 1, 0), 0);
        for (auto& term : terms){
            int d = term.first, c = term.second;
            if (d > 0 && d < n) q[d - 1] = (long long)d * c % mod;
        }
        for (int m = 0; m < n - 1; m++){
            long long sum = q[m];
            for (int j = 1; j < (int)terms.size() && terms[j].first <= m; j++){
                sum = (sum + (long long)(mod - terms[j].second) * q[m - terms[j].first]) % mod;
            }
            q[m] = sum;
        }

        for (int i = 1; i < n; i++) res[i] = q[i - 1] * inv[i] % mod;
        return res;
    }

    vector<int> pow(long long k, int n) const{
        assert(k >= 0 && n >= 0);
        vector<int> res(n);
        if (n == 0) return res;
        if (k == 0){
            res[0] = 1;
            return res;
        }
        if (terms.empty() || terms[0].first >= n) return res;

        int s = terms[0].first, c0 = terms[0].second;
        if (s > 0 && k >= (n + s - 1) / s) return res;

        int shift = s * k;
        vector<int> root = unit_power(k % mod, n - shift, n);
        long long scale = power(c0, k);
        for (int i = shift; i < n; i++) res[i] = root[i - shift] * scale % mod;
        return res;
    }

    /// (f / (c x^s))^e mod x^len for the lowest term c x^s, e a residue mod p, from T g' = e T' g with T(0) = 1
    vector<int> unit_power(long long e, int len, int n) const{
        int s = terms[0].first, top = min(s + len, n);
        long long c0_inv = power(terms[0].second, mod - 2);
        vector<pair<int, long long>> t;
        for (auto& term : terms){
            int d = term.first, c = term.second;
            if (d > s && d < top) t.push_back({d - s, c * c0_inv % mod});
        }

        vector<int> res(len), inv = inverses(len);
        res[0] = 1;
        for (int i = 1; i < len; i++){
            long long sum = 0;
            for (auto& term : t){
                int d = term.first;
                long long c = term.second;
                if (d > i) break;
                long long w = (e * d - i + d) % mod;
                if (w < 0) w += mod;
                sum = (sum + c * w % mod * res[i - d]) % mod;
            }
            res[i] = sum * inv[i] % mod;
        }

        return res;
    }

    /// inv[i] = 1 / i for 1 <= i < n, through inv[i] = -(mod / i) inv[mod % i], which needs a prime mod > i
    vector<int> inverses(int n) const{
        assert(n <= mod);
        vector<int> inv(max(n, 2), 1);
        for (int i = 2; i < n; i++) inv[i] = (long long)(mod - mod / i) * inv[mod % i] % mod;
        return inv;
    }

    long long power(long long b, long long e) const{
        long long res = 1;
        b %= mod;

        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }

        return res;
    }

    /// Tonelli-Shanks, -1 when x is not a quadratic residue, otherwise the smaller of the two roots
    int mod_sqrt(int x) const{
        if (power(x, (mod - 1) / 2) != 1) return -1;

        int s = 0, q = mod - 1, z = 2;
        while (q % 2 == 0){
            q /= 2;
            s++;
        }
        while (power(z, (mod - 1) / 2) != mod - 1) z++;

        long long c = power(z, q), t = power(x, q), r = power(x, (q + 1) / 2);
        while (t != 1){
            int i = 0;
            for (long long u = t; u != 1; u = u * u % mod) i++;

            long long b = c;
            for (int j = 0; j < s - i - 1; j++) b = b * b % mod;
            s = i;
            c = b * b % mod;
            t = t * c % mod;
            r = r * b % mod;
        }

        return min<long long>(r, mod - r);
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
        while (isdigit(c)){
            v = v * 10 - (c - '0');

            /// The digits already buffered are read through a local pointer, 15% faster than a read_char call per digit
            const char *p = inbuf + inptr, *end = inbuf + buf_len;
            while (p < end && (unsigned)(*p - '0') < 10) v = v * 10 - (*p++ - '0');
            inptr = p - inbuf;

            /// A buffered non-digit ends the number and is eaten, else the digits may go on after a refill
            if (p < end){
                inptr++;
                break;
            }
            c = read_char();
        }
        x = neg ? v : -v;

        return true;
    }

    /// The C locale's isspace set, false for EOF, and 1.4-1.7x faster for strings and chars than isspace's table lookup
    inline bool is_space(int c){
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (is_space(c)) c = read_char();
        if (c == EOF) return 0;

        for (; !is_space(c) && c != EOF; c = read_char()){
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

        while (is_space(c)) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !is_space(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (is_space(c));
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

#define MAX 7490

const int coins[] = {1, 5, 10, 25, 50};

int main(){
    int n, i;
    map<int, long long> prod;
    vector<pair<int, long long>> terms;

    prod[0] = 1;
    for (i = 0; i < 5; i++){
        map<int, long long> next = prod;
        for (auto it : prod) next[it.first + coins[i]] -= it.second;
        prod = next;
    }
    for (auto it : prod) terms.push_back(it);

    vector<int> ways = SparseSeries(terms, 2147483647).inverse(MAX);
    while (read(n)) write(ways[n]);

    return 0;
}
