/***
 *
 * Modular Roots
 * Multiplicative order, primitive root, square root and k-th roots modulo a prime
 *
 * Complexity:
 *   multiplicative_order(a, m): O(m^(1/4) log m) expected to factor m and phi(m), then O(log^2 m)
 *   primitive_root(m): the same factoring, then O(log^2 m) per candidate tried
 *   sqrt_mod(a, p): O(log^2 p + c log p) for c the smallest quadratic non-residue mod p, below 2 ln^2 p under GRH
 *   kth_root(a, k, p): O(p^(1/4) log p) expected to factor g = gcd(k, p - 1), then O(log^2 p) per prime of g
 *   kth_roots(a, k, p): kth_root, plus O(g log g) to list and sort the g roots
 *
 * kth_root(a, k, p): any x with x^k = a mod p, -1 if there is none
 * kth_roots(a, k, p): every such x in increasing order, gcd(k, p - 1) of them when a != 0 has a root
 * multiplicative_order(a, m): smallest x >= 1 with a^x = 1 mod m, -1 if gcd(a, m) != 1
 * primitive_root(m): smallest g whose order is phi(m), -1 if none exists
 * sqrt_mod(a, p): the smaller of the two x with x^2 = a mod p (the other is p - x), -1 if there is none
 *
 * Notes:
 *   - p must be prime (not checked), every modulus in [1, 2^62), k >= 0, a is any long long and is reduced mod m
 *   - 0^0 = 1, so kth_roots(1, 0, p) returns all p residues
 *   - A primitive root exists only for m = 1, 2, 4, p^k and 2p^k with p an odd prime, primitive_root(1) = 0
 *   - kth_root adds O(sqrt(q) log p) time and O(sqrt(q)) memory per prime q dividing both g and (p - 1) / g
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Factors with its own small Pollard rho, rho::factorize and rho::is_prime, not the tuned pollard_rho.cpp
 *   - sqrt_mod alone is in modular_sqrt.cpp, a much shorter paste
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY modular_sqrt from code_library/number_theory/modular_sqrt.cpp
long long mul_mod(long long a, long long b, long long m){
    return (unsigned __int128)a * b % m;
}

long long pow_mod(long long x, long long e, long long m){
    long long res = 1 % m;
    for (x %= m; e; e >>= 1, x = mul_mod(x, x, m)){
        if (e & 1) res = mul_mod(res, x, m);
    }
    return res;
}

/// Smallest n that is not a q-th power mod p, for a prime q dividing p - 1
long long non_residue(long long q, long long p){
    long long n = 2;
    while (pow_mod(n, (p - 1) / q, p) == 1) n++;
    return n;
}

long long sqrt_mod(long long a, long long p){
    a = (a % p + p) % p;
    if (a == 0) return 0;
    if (pow_mod(a, (p - 1) / 2, p) != 1) return -1;

    long long s = p - 1;
    int r = 0;
    for (; s % 2 == 0; s /= 2) r++;

    /// Invariant: x^2 = a * b, and the order of b is 2^m for some m < r
    long long x = pow_mod(a, (s + 1) / 2, p), b = pow_mod(a, s, p);
    /// Only b != 1 needs a non-residue: p = 2 has none, and skipping the search is 21-27% faster on random 62-bit p
    long long g = b == 1 ? 1 : pow_mod(non_residue(2, p), s, p);
    while (b != 1){
        int m = 0;
        for (long long t = b; t != 1; t = mul_mod(t, t, p)) m++;

        long long gs = pow_mod(g, 1LL << (r - m - 1), p);
        g = mul_mod(gs, gs, p), x = mul_mod(x, gs, p), b = mul_mod(b, g, p), r = m;
    }

    return min(x, p - x);
}
/// END COPY modular_sqrt

/// Requires gcd(a, m) = 1, returns 0 for m = 1
long long inverse_mod(long long a, long long m){
    long long r = a % m, next_r = m, s = 1, next_s = 0;

    while (next_r){
        long long q = r / next_r;
        r -= q * next_r, swap(r, next_r);
        s -= q * next_s, swap(s, next_s);
    }

    return (s % m + m) % m;
}

namespace rho{
    bool is_prime(long long n){
        static const long long bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

        if (n < 2) return false;
        for (long long q : bases){
            if (n % q == 0) return n == q;
        }

        long long d = n - 1;
        int s = 0;
        for (; !(d & 1); d >>= 1) s++;

        for (long long base : bases){
            long long x = pow_mod(base, d, n);
            if (x == 1) continue;
            for (int i = 1; i < s && x != n - 1; i++) x = mul_mod(x, x, n);
            if (x != n - 1) return false;
        }

        return true;
    }

    /// Brent's cycle detection with one gcd per 128 steps, returns a proper divisor of a composite n
    long long pollard_rho(long long n){
        if (n % 2 == 0) return 2;

        for (long long c = 1; ; c++){
            auto next = [&](long long x){ return (mul_mod(x, x, n) + c) % n; };
            long long x = 0, y = 0, ys = 0, q = 1, g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y);
                for (long long k = 0; k < r && g == 1; k += 128){
                    ys = y;
                    for (long long i = 0; i < 128 && i < r - k; i++){
                        y = next(y);
                        q = mul_mod(q, abs(x - y), n);
                    }
                    g = gcd(q, n);
                }
            }

            /// The batch overshot to a multiple of n, replay it one step at a time
            if (g == n){
                do{
                    ys = next(ys);
                    g = gcd(abs(x - ys), n);
                } while (g == 1);
            }
            if (g != n) return g;
        }
    }

    void collect_primes(long long n, vector<long long>& primes){
        if (n == 1) return;
        if (is_prime(n)){
            primes.push_back(n);
            return;
        }

        long long d = pollard_rho(n);
        collect_primes(d, primes);
        collect_primes(n / d, primes);
    }

    /// Sorted (prime, exponent) pairs
    vector<pair<long long, int>> factorize(long long n){
        vector<long long> primes;

        collect_primes(n, primes);
        sort(primes.begin(), primes.end());

        vector<pair<long long, int>> res;
        for (long long q : primes){
            if (!res.empty() && res.back().first == q) res.back().second++;
            else res.emplace_back(q, 1);
        }

        return res;
    }
}

/// Some x with x^(q^e) = b mod p, for a prime power q^e dividing p - 1 and b != 0 a q^e-th power
long long prime_power_root(long long b, long long q, int e, long long p){
    long long t = p - 1, qe = 1;
    int s = 0;

    for (; t % q == 0; t /= q) s++;
    for (int i = 0; i < e; i++) qe *= q;

    /// x^(q^e) / b = (b^t)^((q^e - 1) / t mod q^e) lies in the subgroup of order q^s, which z generates
    long long d = ((unsigned __int128)(qe - 1) * inverse_mod(t % qe, qe) % qe * t + 1) / qe;
    long long x = pow_mod(b, d, p), err = mul_mod(b, inverse_mod(pow_mod(x, qe, p), p), p);
    long long z = pow_mod(non_residue(q, p), t, p), w = pow_mod(z, (p - 1) / t / q, p);

    long long m = sqrtl(q) + 1, giant = 0;
    vector<pair<long long, long long>> baby;

    /// Baby step giant step in the order q subgroup of w, built lazily as h is always 1 when q^e is all of q in p - 1
    auto log_w = [&](long long h){
        if (h == 1) return 0LL;
        if (baby.empty()){
            baby.resize(m);
            for (long long j = 0, cur = 1; j < m; j++, cur = mul_mod(cur, w, p)) baby[j] = {cur, j};
            sort(baby.begin(), baby.end());
            giant = inverse_mod(pow_mod(w, m, p), p);
        }

        for (long long i = 0; ; i++, h = mul_mod(h, giant, p)){
            auto it = lower_bound(baby.begin(), baby.end(), make_pair(h, 0LL));
            if (it != baby.end() && it->first == h) return i * m + it->second;
        }
    };

    /// Pohlig-Hellman: err = z^log_err, one base q digit per step, the lowest e digits come out 0
    long long z_inv = inverse_mod(z, p), log_err = 0;
    for (long long i = 0, qi = 1; i < s; i++, qi *= q){
        long long h = mul_mod(err, pow_mod(z_inv, log_err, p), p);
        log_err += log_w(pow_mod(h, (p - 1) / t / qi / q, p)) * qi;
    }

    return mul_mod(x, pow_mod(z, log_err / qe, p), p);
}

long long kth_root(long long a, long long k, long long p){
    a = (a % p + p) % p;
    if (k == 0) return a == 1 ? 0 : -1;
    if (a == 0) return 0;

    long long g = gcd(k, p - 1), m = (p - 1) / g;
    if (pow_mod(a, m, p) != 1) return -1;

    /// a is a g-th power, so b = a^(1 / (k / g)) is one too, and any g-th root of b is a k-th root of a
    long long b = pow_mod(a, inverse_mod(k / g % m, m), p);
    for (auto [q, e] : rho::factorize(g)) b = prime_power_root(b, q, e, p);
    return b;
}

vector<long long> kth_roots(long long a, long long k, long long p){
    a = (a % p + p) % p;
    if (k == 0 && a == 1){
        vector<long long> all(p);
        iota(all.begin(), all.end(), 0);
        return all;
    }

    long long x = kth_root(a, k, p);
    if (x <= 0) return x ? vector<long long>{} : vector<long long>{0};

    /// w has order exactly g, so x * w^i for i < g are the g distinct roots
    long long g = gcd(k, p - 1), w = 1;
    for (auto [q, e] : rho::factorize(g)){
        long long cofactor = p - 1;
        for (int i = 0; i < e; i++) cofactor /= q;
        w = mul_mod(w, pow_mod(non_residue(q, p), cofactor, p), p);
    }

    vector<long long> roots(g);
    for (long long i = 0; i < g; i++, x = mul_mod(x, w, p)) roots[i] = x;
    sort(roots.begin(), roots.end());
    return roots;
}

long long multiplicative_order(long long a, long long m){
    a = (a % m + m) % m;
    if (gcd(a, m) != 1) return -1;

    long long phi = m;
    for (auto [q, e] : rho::factorize(m)) phi = phi / q * (q - 1);

    long long order = phi;
    for (auto [q, e] : rho::factorize(phi)){
        while (order % q == 0 && pow_mod(a, order / q, m) == 1) order /= q;
    }

    return order;
}

long long primitive_root(long long m){
    if (m <= 4) return m - 1;

    long long odd = m % 2 ? m : m / 2;
    auto f = rho::factorize(odd);
    if (m % 4 == 0 || f.size() != 1) return -1;

    long long phi = odd / f[0].first * (f[0].first - 1);
    auto phi_factors = rho::factorize(phi);
    for (long long g = 2; ; g++){
        if (gcd(g, m) != 1) continue;
        bool is_root = true;
        for (auto [q, e] : phi_factors) is_root = is_root && pow_mod(g, phi / q, m) != 1;
        if (is_root) return g;
    }
}

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

#define MAX 105
#define clr(ar) memset(ar, 0, sizeof(ar))

const long long P = 998244353;
const long long MOD = P - 1;

int k;
long long b[MAX], T[MAX][MAX], R[MAX][MAX], tmp[MAX][MAX];

void multiply(long long A[MAX][MAX], long long B[MAX][MAX]){
    int i, j, l;

    clr(tmp);
    for (i = 0; i < k; i++){
        for (l = 0; l < k; l++){
            for (j = 0; j < k; j++) tmp[i][j] = (tmp[i][j] + A[i][l] * B[l][j]) % MOD;
        }
    }
    memcpy(A, tmp, sizeof(tmp));
}

long long exponent(long long n){
    int i;

    clr(T), clr(R);
    for (i = 0; i < k; i++) R[i][i] = 1, T[0][i] = b[i];
    for (i = 1; i < k; i++) T[i][i - 1] = 1;

    for (; n; n >>= 1){
        if (n & 1) multiply(R, T);
        multiply(T, T);
    }

    return R[0][0];
}

int main(){
    int i;
    long long n, m, e;

    read(k);
    for (i = 0; i < k; i++) read(b[i]);
    read(n, m);

    e = exponent(n - k);
    if (!e) write(m == 1 ? 1 : -1);
    else write(kth_root(m, e, P));
    return 0;
}
