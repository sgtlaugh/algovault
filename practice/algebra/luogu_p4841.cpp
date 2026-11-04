/***
 *
 * Number Theoretic Transform
 * Polynomial multiplication modulo an NTT-friendly prime, modulo any m < 2^31, or exact in 64 bits
 *
 * Complexity: O(n log n) for a product of total length n
 *
 * ntt::crt_multiply(a, b): a * b modulo M = 167772161 * 469762049 * 754974721 ~ 5.95e25, unsigned __int128 in [0, M)
 * ntt::exact_multiply(a, b): the exact a * b of any long long inputs while every true coefficient is in [-2^63, 2^63)
 * ntt::mod_multiply(a, b, m): a * b modulo any 1 <= m < 2^31, prime or not
 * ntt::multiply<MOD, ROOT>(a, b): a * b modulo MOD, default MOD = 998244353, ROOT found from MOD
 *
 * Notes:
 *   - The result length is |a| + |b| - 1, at most 2^24 for the other three, and an empty input gives an empty result
 *   - multiply needs the result length rounded up to a power of two to divide MOD - 1
 *   - multiply allows result lengths up to 2^23 modulo 998244353 and 2^20 modulo 7340033 = 7 * 2^20 + 1
 *   - MOD must be a prime below 2^31, an explicit ROOT a quadratic non-residue modulo MOD, which any primitive root is
 *   - Inputs may be negative or >= the modulus, multiply and mod_multiply return values in [0, modulus)
 *   - Roots stay cached per modulus, 8 bytes per element of the largest padded length, three moduli for the others
 *   - Requires __int128 (64-bit GCC or Clang) for exact_multiply and crt_multiply
 *   - Use fft.cpp for exact products with max(n, m) * max|value|^2 below ~4.5e13: one transform there, three here
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace ntt{
    /// BEGIN SHARED ntt_multiply
    const unsigned P1 = 167772161, P2 = 469762049, P3 = 754974721;

    constexpr bool is_prime(unsigned n){
        if (n < 2) return false;
        for (unsigned long long d = 2; d * d <= n; d++){
            if (n % d == 0) return false;
        }

        return true;
    }

    /// BEGIN SHARED ntt_core
    constexpr unsigned long long power(unsigned long long b, unsigned long long e, unsigned long long mod){
        unsigned long long res = 1;
        b %= mod;
        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }

        return res;
    }

    /// Newton's step y (2 - p y) doubles the correct low bits of 1 / p modulo 2^32, y = p starts with 3 of them
    constexpr unsigned inverse_mod_2_32(unsigned p, unsigned y, int steps){
        return steps ? inverse_mod_2_32(p, y * (2 - p * y), steps - 1) : y;
    }

    /// Values live in [0, BOUND), 2 MOD below 2^30 to skip the final reduction, else MOD, so u - v + BOUND < 2^32
    template<unsigned MOD>
    constexpr unsigned bound(){
        return MOD < (1u << 30) ? 2 * MOD : MOD;
    }

    /// Montgomery reduction, x / 2^32 modulo MOD for x < MOD 2^32, in [0, bound<MOD>()), replacing the slower x % MOD
    template<unsigned MOD>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = 0u - inverse_mod_2_32(MOD, MOD, 4);
        unsigned r = (x + (unsigned long long)((unsigned)x * NEG_INV) * MOD) >> 32;
        return MOD < (1u << 30) || r < MOD ? r : r - MOD;
    }

    /// rt[k + j] = w^j R modulo MOD for each power of two k < n, w a primitive 2k-th root of unity or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), ((unsigned long long)1 << 32) % MOD);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            z = (z << 32) % MOD;
            for (int i = k; i < 2 * k; i++){
                unsigned w = redc<MOD>(rt[i >> 1] * z);
                rt[i] = i & 1 ? (w >= MOD ? w - MOD : w) : rt[i >> 1];
            }
        }

        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two at most rt.size()
    template<unsigned MOD>
    void forward_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = redc<MOD>((unsigned long long)(u - v + BOUND) * rt[j + k]);
                }
            }
        }
    }

    /// Decimation in time with the inverse roots, bit-reversed order in, natural order out, scaled by n
    template<unsigned MOD>
    void inverse_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<MOD>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = u >= v ? u - v : u + BOUND - v;
                }
            }
        }
    }
    /// END SHARED ntt_core

    /// Any quadratic non-residue g works as ROOT: g^((MOD - 1) / n) has order n for each power of two n | MOD - 1
    constexpr unsigned non_residue(unsigned mod){
        unsigned g = 2;
        while (power(g, (mod - 1) / 2, mod) != mod - 1) g++;
        return g;
    }

    template<unsigned MOD = 998244353, unsigned ROOT = non_residue(MOD)>
    vector<long long> multiply(const vector<long long>& a, const vector<long long>& b){
        static_assert(MOD < (1u << 31) && is_prime(MOD), "MOD must be a prime below 2^31");
        static_assert(power(ROOT, (MOD - 1) / 2, MOD) == MOD - 1, "ROOT must be a quadratic non-residue modulo MOD");
        if (a.empty() || b.empty()) return {};

        int s = a.size() + b.size() - 1, n = 1;
        while (n < s) n <<= 1;
        assert((MOD - 1) % n == 0);

        /// Each level's roots sit at fixed positions, so one table grown to the largest n serves every smaller call
        static vector<unsigned> rt, irt;
        if ((int)rt.size() < n) rt = roots<MOD, ROOT>(n, false), irt = roots<MOD, ROOT>(n, true);

        bool same = a == b;
        vector<unsigned> x(n, 0), y;
        for (size_t i = 0; i < a.size(); i++) x[i] = (a[i] % (long long)MOD + MOD) % MOD;
        forward_transform<MOD>(x, rt);
        if (!same){
            y.assign(n, 0);
            for (size_t i = 0; i < b.size(); i++) y[i] = (b[i] % (long long)MOD + MOD) % MOD;
            forward_transform<MOD>(y, rt);
        }

        const vector<unsigned>& z = same ? x : y;
        unsigned long long inv = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * z[i] % MOD * inv % MOD;
        inverse_transform<MOD>(x, irt);

        vector<long long> res(s);
        for (int i = 0; i < s; i++) res[i] = x[i] >= MOD ? x[i] - MOD : x[i];
        return res;
    }
    /// END SHARED ntt_multiply

    vector<unsigned __int128> crt_multiply(const vector<long long>& a, const vector<long long>& b){
        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);

        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        vector<unsigned __int128> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            /// Garner: res = x1 + x2 * P1 + x3 * P1 * P2, each digit chosen so res matches c1, c2, c3 modulo P1, P2, P3
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = x1 + (unsigned __int128)x2 * P1 + (unsigned __int128)x3 * P1 * P2;
        }

        return res;
    }

    vector<long long> exact_multiply(const vector<long long>& a, const vector<long long>& b){
        const unsigned __int128 M = (unsigned __int128)P1 * P2 * P3;

        auto c = crt_multiply(a, b);
        vector<long long> res(c.size());

        /// A negative true value v comes back as M + v, which lies above M / 2
        for (size_t i = 0; i < c.size(); i++) res[i] = c[i] > M / 2 ? (long long)((__int128)c[i] - (__int128)M) : (long long)c[i];
        return res;
    }

    /// BEGIN SHARED ntt_mod_multiply
    vector<long long> mod_multiply(vector<long long> a, vector<long long> b, long long m){
        assert(1 <= m && m < (1LL << 31));
        for (auto& x : a) x = (x % m + m) % m;
        for (auto& x : b) x = (x % m + m) % m;

        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (unsigned long long)P1 * P2 % m;
        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        /// True coefficients < 2^23 m^2 < M ~ 5.95e25 make x1 + x2 P1 + x3 P1 P2 exact, reduced below 2^62 per term
        vector<long long> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
        }

        return res;
    }
    /// END SHARED ntt_mod_multiply
}

/***
 *
 * Online Convolution
 * A sequence f whose every term depends on the convolution of the earlier terms with a sequence g known in advance
 *
 * Complexity: O(n log^2 n) with an O(n log n) multiply, O(n) calls to it
 *
 * online_convolution(n, g, mod, multiply, next): f_i = next(i, s_i) for i < n, s_i = sum over j < i of f_j g_(i-j)
 *
 * Notes:
 *   - multiply(a, b) returns the product of two vector<long long> modulo mod in [0, mod), like ntt::multiply
 *   - multiply may drop trailing zeros, like Polynomial's operator*, so ntt::mod_multiply or an FFT plug in as well
 *   - next is called once per i in increasing order with s_i in [0, mod), so it may keep state between calls
 *   - 1 <= mod < 2^31, g_0 is never read, missing g terms are 0, and values of g and next are reduced into [0, mod)
 *   - g = f itself, as in f_i = next(i, (f * f)_i), is not supported: g must be known before f is built
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED online_convolution
template<class Multiply, class Next>
vector<long long> online_convolution(int n, vector<long long> g, long long mod, Multiply multiply, Next next){
    assert(n >= 0 && 1 <= mod && mod < (1LL << 31));
    const int LEAF = 32;

    g.resize(n, 0);
    for (auto& x : g) x = (x % mod + mod) % mod;
    vector<long long> f(n, 0), s(n, 0);

    auto solve = [&](auto& self, int l, int r) -> void{
        if (r - l <= LEAF){
            for (int i = l; i < r; i++){
                for (int j = l; j < i; j++) s[i] = (s[i] + f[j] * g[i - j]) % mod;
                f[i] = (next(i, s[i]) % mod + mod) % mod;
            }

            return;
        }

        int mid = (l + r) / 2;
        self(self, l, mid);

        /// f_j for j in [l, mid) times g_(i-j) lands on index i - l of the product, for every i in [mid, r)
        vector<long long> prod = multiply(vector<long long>(f.begin() + l, f.begin() + mid),
                                          vector<long long>(g.begin(), g.begin() + (r - l)));
        /// A multiply that trims trailing zeros, like Polynomial's, returns fewer terms
        for (int i = mid; i < r && i - l < (int)prod.size(); i++) s[i] = (s[i] + prod[i - l]) % mod;

        self(self, mid, r);
    };

    solve(solve, 0, n);
    return f;
}
/// END SHARED online_convolution

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

const int MOD = 1004535809;

#define MAX 130010

long long fact[MAX], inv_fact[MAX], graphs[MAX];

long long expo(long long b, long long e){
    long long res = 1;
    for (b %= MOD; e; e >>= 1, b = b * b % MOD){
        if (e & 1) res = res * b % MOD;
    }

    return res;
}

long long next_term(int i, long long s){
    return i ? (graphs[i] * fact[i] % MOD * inv_fact[i - 1] - s) % MOD : 0;
}

int main(){
    int n, i;

    read(n);
    fact[0] = 1;
    for (i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % MOD;
    inv_fact[n] = expo(fact[n], MOD - 2);
    for (i = n; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % MOD;

    vector<long long> g(n + 1);
    for (i = 0; i <= n; i++){
        graphs[i] = expo(2, (long long)i * (i - 1) / 2) * inv_fact[i] % MOD;
        g[i] = graphs[i];
    }

    vector<long long> c = online_convolution(n + 1, g, MOD, ntt::multiply<MOD>, next_term);
    write(c[n] * fact[n - 1] % MOD);
    return 0;
}
