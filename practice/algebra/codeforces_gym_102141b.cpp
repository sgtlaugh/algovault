/***
 *
 * Pollard Rho
 * Integer factorization and deterministic primality test for 1 <= n < 2^63
 *
 * Complexity: expected O(n^(1/4)) per split, O(1e6) once to build the sieve
 *
 * PollardRho f;           sieves smallest prime factors up to 1e6, keep a single instance around
 * f.factorize(n):         prime factors of n in ascending order with multiplicity, empty for n = 1
 * f.is_prime(n):          deterministic Miller Rabin for n < 2^63
 *
 * Brent's cycle detection with one gcd per batch of 128 steps and Montgomery multiplication
 * Requires __int128 (64-bit GCC or Clang)
 * pisano_period.cpp carries checked copies of Montgomery and PollardRho, miller_rabin.cpp of Montgomery
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n
/// BEGIN SHARED montgomery
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// Newton iteration, doubles the correct bits of n^-1 mod 2^64 each step
        r2 = -n % n;
        r2 = (unsigned __int128)r2 * r2 % n;
    }

    unsigned long long reduce(unsigned __int128 x) const{
        unsigned long long q = (unsigned long long)x * inv, m = ((unsigned __int128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    unsigned long long mul(unsigned long long x, unsigned long long y) const{
        return reduce((unsigned __int128)x * y);
    }

    unsigned long long add(unsigned long long x, unsigned long long y) const{
        return (x += y) >= n ? x - n : x;
    }

    unsigned long long to_mont(unsigned long long x) const{
        return mul(x, r2);
    }

    unsigned long long pow(unsigned long long x, unsigned long long e) const{
        unsigned long long res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }
        return res;
    }
};
/// END SHARED montgomery

/// BEGIN SHARED pollard_rho
struct PollardRho{
    static constexpr int LIMIT = 1000001;
    vector<int> spf;
    mt19937_64 rng;

    PollardRho() : spf(LIMIT, 0), rng(20260105){
        vector<int> primes;
        for (int i = 2; i < LIMIT; i++){
            if (!spf[i]) spf[i] = i, primes.push_back(i);
            for (int p : primes){
                if (p > spf[i] || (long long)i * p >= LIMIT) break;
                spf[i * p] = p;
            }
        }
    }

    bool is_prime(long long n) const{
        if (n < LIMIT) return n >= 2 && spf[n] == n;
        if (!(n & 1)) return false;

        static const long long BASES[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
        long long s = 0, r = n - 1;
        for (; !(r & 1); r >>= 1, s++) {}

        Montgomery mont(n);
        unsigned long long one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        for (long long base : BASES){
            long long a = base % n;
            if (!a) continue;  /// n divides the base, this base says nothing about n
            unsigned long long c = mont.pow(mont.to_mont(a), r);
            for (int j = 0; j < s; j++){
                unsigned long long d = mont.mul(c, c);
                if (d == one && c != one && c != minus_one) return false;
                c = d;
            }

            if (c != one) return false;
        }
        return true;
    }

    vector<long long> factorize(long long n){
        assert(n >= 1);
        vector<long long> res;
        collect(n, res);
        sort(res.begin(), res.end());
        return res;
    }

    void collect(long long n, vector<long long>& out){
        if (n < LIMIT){
            for (; n > 1; n /= spf[n]) out.push_back(spf[n]);
            return;
        }

        if (!(n & 1)){
            out.push_back(2);
            return collect(n >> 1, out);
        }

        if (is_prime(n)){
            out.push_back(n);
            return;
        }

        long long d = split(n);
        collect(d, out);
        collect(n / d, out);
    }

    /// A nontrivial factor of an odd composite n
    long long split(long long n){
        const int BATCH = 128;
        Montgomery mont(n);
        auto next = [&](unsigned long long x, unsigned long long c){ return mont.add(mont.mul(x, x), c); };
        auto dist = [](unsigned long long x, unsigned long long y){ return x > y ? x - y : y - x; };

        while (true){
            unsigned long long c = rng() % (n - 1) + 1, y = rng() % n, x = y, ys = y, q = mont.to_mont(1), g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y, c);

                for (long long k = 0; k < r && g == 1; k += BATCH){
                    ys = y;
                    for (long long i = 0; i < BATCH && i < r - k; i++){
                        y = next(y, c);
                        q = mont.mul(q, dist(x, y));
                    }
                    g = __gcd(q, (unsigned long long)n);
                }
            }

            if (g == (unsigned long long)n){  /// the batch overshot to a multiple of n, replay it one step at a time
                do{
                    ys = next(ys, c);
                    g = __gcd(dist(x, ys), (unsigned long long)n);
                } while (g == 1);
            }

            if (g != (unsigned long long)n) return g;
        }
    }
};
/// END SHARED pollard_rho

/***
 *
 * Walsh Hadamard Transform
 * xor, or and and convolutions of two vectors, exact or modulo mod
 *
 * Complexity: O(n log n)
 *
 * WalshHadamard fwht: a convolver for any number of calls of any length
 * fwht.xor_convolution(A, B, mod = 0): C where C[k] is the sum of A[i] * B[j] over all i ^ j = k, reduced modulo mod when mod > 0
 * fwht.or_convolution(A, B, mod = 0): the same with i | j = k
 * fwht.and_convolution(A, B, mod = 0): the same with i & j = k
 * WalshHadamard::walsh_transform(ar, n, conv_type, mod = 0): transforms ar[0..n) in place, conv_type WalshHadamard::OR, AND or XOR
 * WalshHadamard::inverse_walsh_transform(ar, n, conv_type, mod = 0): undoes walsh_transform
 * With a mod, the two transforms need every ar[i] already in [0, mod)
 *
 * A and B must have the same length n, and n must be a power of 2
 * Exact mode (mod = 0): no overflow while n^2 * max|A[i]| * max|B[j]| < 2^62
 * Modular mode: 1 <= mod < 2^31, any input values (negative or >= mod are reduced), results in [0, mod)
 * xor needs an odd mod, or and and take any mod
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct WalshHadamard{
    enum ConvType{OR, AND, XOR};

    vector<long long> buffer;  /// reused across calls: allocating a fresh 8 MB vector per call at n = 2^20 costs ~10% in page faults

    vector<long long> xor_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, XOR, mod);
    }

    vector<long long> or_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, OR, mod);
    }

    vector<long long> and_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, AND, mod);
    }

    static void walsh_transform(long long* ar, int n, ConvType conv_type, long long mod = 0){
        dispatch<false>(ar, n, conv_type, mod);
    }

    static void inverse_walsh_transform(long long* ar, int n, ConvType conv_type, long long mod = 0){
        dispatch<true>(ar, n, conv_type, mod);
    }

    vector<long long> convolution(const vector<long long>& A, const vector<long long>& B, ConvType conv_type, long long mod = 0){
        int n = A.size();
        assert(A.size() == B.size() && __builtin_popcount(n) == 1);
        vector<long long> res = A;
        buffer.assign(B.begin(), B.end());
        if (mod){
            for (auto& x : res) if ((x %= mod) < 0) x += mod;
            for (auto& x : buffer) if ((x %= mod) < 0) x += mod;
        }

        walsh_transform(res.data(), n, conv_type, mod);
        walsh_transform(buffer.data(), n, conv_type, mod);
        if (mod) for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i] % mod;
        else for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i];
        inverse_walsh_transform(res.data(), n, conv_type, mod);

        return res;
    }

    /// conv_type and mod != 0 become template arguments so the butterfly loops carry no per-element branches
    template<bool inverse>
    static void dispatch(long long* ar, int n, ConvType conv_type, long long mod){
        assert(mod >= 0 && mod < (1LL << 31) && (conv_type != XOR || !mod || mod % 2 == 1));
        if (conv_type == OR) mod ? transform<OR, inverse, true>(ar, n, mod) : transform<OR, inverse, false>(ar, n, mod);
        if (conv_type == AND) mod ? transform<AND, inverse, true>(ar, n, mod) : transform<AND, inverse, false>(ar, n, mod);
        if (conv_type == XOR) mod ? transform<XOR, inverse, true>(ar, n, mod) : transform<XOR, inverse, false>(ar, n, mod);
    }

    template<ConvType conv_type, bool inverse, bool modular>
    static void transform(long long* ar, int n, long long mod){
        if (n < 2) return;

        int m = n >> 1;
        transform<conv_type, inverse, modular>(ar, m, mod);
        transform<conv_type, inverse, modular>(ar + m, m, mod);

        for (int i = 0; i < m; i++){
            long long x = ar[i], y = ar[i + m];
            if constexpr (conv_type == OR) ar[i + m] = inverse ? sub<modular>(y, x, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == AND) ar[i] = inverse ? sub<modular>(x, y, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == XOR){
                ar[i] = add<modular>(x, y, mod), ar[i + m] = sub<modular>(x, y, mod);
                if constexpr (inverse) ar[i] = half<modular>(ar[i], mod), ar[i + m] = half<modular>(ar[i + m], mod);
            }
        }
    }

    /// modular operands are in [0, mod); (s >> 63) & mod adds mod back exactly when s went negative
    template<bool modular>
    static long long add(long long x, long long y, long long mod){
        long long s = x + y;
        if constexpr (modular) s -= mod, s += (s >> 63) & mod;
        return s;
    }

    template<bool modular>
    static long long sub(long long x, long long y, long long mod){
        long long d = x - y;
        if constexpr (modular) d += (d >> 63) & mod;
        return d;
    }

    /// exact mode: v is even, so >> 1 divides; modular mode: mod is odd, so adding it to an odd v makes v + mod even
    template<bool modular>
    static long long half(long long v, long long mod){
        if constexpr (modular) v += (v & 1) * mod;
        return v >> 1;
    }
};

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

#define MOD 2147483647

int main(){
    int t, T = 0, n, i, j, k, mask, res;
    long long X, a;
    char buf[64];
    PollardRho rho;

    read(t);
    while (t--){
        read(n, X);
        vector<long long> primes = rho.factorize(X), powers;
        for (i = 0; i < (int)primes.size(); i++){
            if (i && primes[i] == primes[i - 1]) powers.back() *= primes[i];
            else powers.push_back(primes[i]);
        }

        k = powers.size();
        vector<long long> cnt(1 << k, 0);
        for (i = 0; i < n; i++){
            read(a);
            if (X % a) continue;
            mask = 0;
            for (j = 0; j < k; j++){
                if (a % powers[j] == 0) mask |= (1 << j);
            }
            cnt[mask]++;
        }

        WalshHadamard::walsh_transform(cnt.data(), 1 << k, WalshHadamard::OR);
        vector<long long> ways(1 << k, 1);

        res = -1;
        for (j = 1; j <= max(k, 1) && res == -1; j++){
            for (i = 0; i < (1 << k); i++) ways[i] = ways[i] * cnt[i] % MOD;
            vector<long long> exact = ways;
            WalshHadamard::inverse_walsh_transform(exact.data(), 1 << k, WalshHadamard::OR, MOD);
            if (exact[(1 << k) - 1]) res = j;
        }

        snprintf(buf, sizeof(buf), "Case %d: %d", ++T, res);
        write(string(buf));
    }
    return 0;
}
