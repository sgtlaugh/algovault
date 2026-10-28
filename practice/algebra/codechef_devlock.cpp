/***
 *
 * Number Theoretic Transform
 * Polynomial multiplication modulo an NTT-friendly prime, modulo any m < 2^31, or exact in 64 bits
 *
 * Complexity: O(n log n) for a product of total length n
 *
 * ntt::multiply<MOD, ROOT>(a, b): a * b modulo MOD, default MOD = 998244353, ROOT found from MOD, the result length rounded up to a power of two must divide MOD - 1
 * ntt::mod_multiply(a, b, m): a * b modulo any 1 <= m < 2^31, prime or not, result length at most 2^24
 * ntt::exact_multiply(a, b): the exact a * b of any long long inputs, correct whenever every true coefficient lies in [-2^63, 2^63), result length at most 2^24
 * ntt::crt_multiply(a, b): every true coefficient of a * b modulo M = 167772161 * 469762049 * 754974721 ~ 5.95e25 as an unsigned __int128 in [0, M), result length at most 2^24
 *
 * The result length is |a| + |b| - 1, and an empty input gives an empty result
 * MOD must be a prime below 2^31; an explicit ROOT must be a quadratic non-residue modulo MOD, which any primitive root is
 * multiply allows result lengths up to 2^23 modulo 998244353 and 2^20 modulo 7340033 = 7 * 2^20 + 1
 * Inputs may be negative or >= the modulus, multiply and mod_multiply return values in [0, modulus)
 * Requires __int128 (64-bit GCC or Clang) for exact_multiply and crt_multiply
 * Use fft.cpp for exact products of small values like 0..100, about 3x faster there; any other mod runs at a similar speed in both
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace ntt{
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

    /// Values live in [0, BOUND): 2 MOD below 2^30 saves the final reduction, MOD above, so u - v + BOUND stays below 2^32
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

    /// rt[k + j] = w^j R modulo MOD for every power of two k < n, where w is a primitive 2k-th root of unity, or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), 1);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % MOD : rt[i >> 1];
        }

        for (auto& w : rt) w = ((unsigned long long)w << 32) % MOD;
        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two no larger than rt.size()
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

    /// Any quadratic non-residue g works as ROOT: g^((MOD - 1) / n) then has order exactly n for every power of two n dividing MOD - 1
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

        vector<unsigned> x(n, 0), y(n, 0);
        for (size_t i = 0; i < a.size(); i++) x[i] = (a[i] % (long long)MOD + MOD) % MOD;
        for (size_t i = 0; i < b.size(); i++) y[i] = (b[i] % (long long)MOD + MOD) % MOD;

        auto rt = roots<MOD, ROOT>(n, false);
        forward_transform<MOD>(x, rt);
        forward_transform<MOD>(y, rt);

        unsigned long long inv = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * y[i] % MOD * inv % MOD;
        inverse_transform<MOD>(x, roots<MOD, ROOT>(n, true));

        vector<long long> res(s);
        for (int i = 0; i < s; i++) res[i] = x[i] >= MOD ? x[i] - MOD : x[i];
        return res;
    }

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

    vector<long long> mod_multiply(vector<long long> a, vector<long long> b, long long m){
        assert(1 <= m && m < (1LL << 31));
        for (auto& x : a) x = (x % m + m) % m;
        for (auto& x : b) x = (x % m + m) % m;

        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (unsigned long long)P1 * P2 % m;
        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        /// Every true coefficient is below 2^23 m^2 < M ~ 5.95e25, so x1 + x2 P1 + x3 P1 P2 is exact, reduced term by term below 2^62
        vector<long long> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
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
}

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
 *
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ... and returns how many it read before the input ended
 * read accepts integers of any width including __int128, float, double, long double, strings, chars and vectors of these
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... separated by single spaces and ends the line, a vector argument on its own line
 * write accepts integers, strings, chars, vectors of these and fio::fixed values
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 * fio::flush(): writes the pending output now, which also happens at normal exit
 *
 * Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 * A char reads one non-space character
 * A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 * Writing a float without fixed() does not compile, so its precision is always chosen
 * Floats print at printf's speed, since the formatting itself dominates; a floating token keeps its first 127 characters
 * Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: stdin and stdout are one stream each per process, so a second reader instance
/// would split the input between two buffers and a second writer would reorder the output
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

    /// Integers: takes every type the plain overloads below miss, so __int128 needs no is_integral specialization, which breaks C++20
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

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
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

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
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

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

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
}

using namespace fio;

#define MOD 998244353

int n, p, mm, w;
long long cnt[55], inv[15010];

vector<long long> dial_power(long long c, int r){
    vector<long long> a(mm + 1, 0), res((mm + 1) * w, 0);

    a[0] = 1, c %= MOD;
    for (int k = 1; k <= mm; k++){
        long long s = 0;
        for (int j = 1; j <= 9 && j <= k; j++) s = (s + ((c * j - k + j) % MOD + MOD) % MOD * a[k - j]) % MOD;
        a[k] = s * inv[k] % MOD;
    }

    for (int k = 0; k <= mm; k++) res[k * w + k * r % p] = a[k];
    return res;
}

void count_weights(){
    int i, j, x = 1 % p;
    vector<int> first(p, -1), seq;

    for (i = 0; i < n && first[x] == -1; i++){
        first[x] = i, cnt[x]++;
        seq.push_back(x);
        x = x * 10 % p;
    }

    if (i < n){
        int s = first[x], len = i - s, full = (n - i) / len, rem = (n - i) % len;
        for (j = s; j < i; j++) cnt[seq[j]] += full + ((j - s) < rem);
    }
}

int main(){
    int i, k, r;

    read(n, p, mm);
    w = 2 * p - 1;
    inv[1] = 1;
    for (i = 2; i <= mm; i++) inv[i] = (MOD - (MOD / i) * inv[MOD % i] % MOD) % MOD;
    count_weights();

    vector<long long> cur((mm + 1) * w, 0);
    cur[0] = 1;
    for (r = 0; r < p; r++){
        if (!cnt[r]) continue;
        vector<long long> res = ntt::multiply(cur, dial_power(cnt[r], r));

        fill(cur.begin(), cur.end(), 0);
        for (i = 0; i < (int)res.size(); i++){
            k = i / w;
            if (k <= mm) cur[k * w + (i % w) % p] = (cur[k * w + (i % w) % p] + res[i]) % MOD;
        }
    }

    vector<long long> ans(mm + 1);
    long long sum = 0;
    for (k = 0; k <= mm; k++){
        sum = (sum + cur[k * w]) % MOD;
        ans[k] = sum;
    }

    write(ans);
    return 0;
}
