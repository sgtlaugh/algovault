/***
 *
 * Bell Numbers
 * B_n counts the splits of n labelled items into non-empty groups, here modulo m for all n up to a bound or one huge n
 *
 * Complexity:
 *   - bell_number: O(p log p log n)
 *   - bell_numbers: O(n^2) time, O(n) memory
 *   - bell_numbers_ntt: O(n log^2 n)
 *
 * bell_number(n, p): B_n modulo a prime p <= 2^20, for 0 <= n <= 1e18
 * bell_numbers(n, m): B_0, ..., B_n modulo any 1 <= m <= 2^62, prime or not
 * bell_numbers_ntt(n, mod = 998244353): B_0, ..., B_n modulo a prime n < mod <= 2^30, n <= 2^22
 *
 * Notes:
 *   - bell_numbers_ntt is about 3x faster modulo 998244353, other primes take three NTTs per product
 *   - bell_number at n = 1e18 takes about 0.4 s for p = 1e4, 4 s for p = 1e5
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace bell_ntt{
    /// BEGIN COPY stirling_ntt from code_library/combinatorics/stirling_numbers.cpp
    const unsigned MOD = 998244353, P1 = 167772161, P2 = 469762049, P3 = 754974721;

    long long power(long long b, long long e, long long mod){
        long long res = 1 % mod;
        for (b %= mod; e; e >>= 1, b = b * b % mod){
            if (e & 1) res = res * b % mod;
        }

        return res;
    }

    /// Newton's step y (2 - p y) doubles the correct low bits of 1 / p modulo 2^32, y = p starts with 3 of them
    constexpr unsigned inverse_mod_2_32(unsigned p, unsigned y, int steps){
        return steps ? inverse_mod_2_32(p, y * (2 - p * y), steps - 1) : y;
    }

    /// Montgomery reduction with R = 2^32, x / R modulo P in [0, 2P) for any x < P R, replacing the slower x % P
    template<unsigned P>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = 0u - inverse_mod_2_32(P, P, 4);
        return (x + (unsigned long long)((unsigned)x * NEG_INV) * P) >> 32;
    }

    /// rt[k + j] = w^j R mod P for every power of two k < n, where w is a primitive 2k-th root of unity, or its inverse
    template<unsigned P, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), 1);

        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (P - 1) / (2 * k), P);
            if (invert) z = power(z, P - 2, P);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % P : rt[i >> 1];
        }

        for (auto& w : rt) w = ((unsigned long long)w << 32) % P;
        return rt;
    }

    /// Bit-reversed order out, values in [0, 2P), P < 2^30 keeps u - v + 2P below 2^32 and redc inputs below P R
    template<unsigned P>
    void forward(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();

        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= 2 * P ? u + v - 2 * P : u + v;
                    a[i + j + k] = redc<P>((unsigned long long)(u - v + 2 * P) * rt[j + k]);
                }
            }
        }
    }

    /// Inverse roots, bit-reversed order in, natural order out scaled by n, values stay in [0, 2P)
    template<unsigned P>
    void inverse(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();

        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<P>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= 2 * P ? u + v - 2 * P : u + v;
                    a[i + j + k] = u + 2 * P - v >= 2 * P ? u - v : u + 2 * P - v;
                }
            }
        }
    }

    /// Cyclic convolution of length n modulo P, inputs longer than n fold onto index i mod n
    template<unsigned P, unsigned ROOT>
    vector<unsigned> cyclic(const vector<long long>& a, const vector<long long>& b, int n){
        vector<unsigned> x(n, 0), y(n, 0);

        for (size_t i = 0; i < a.size(); i++) x[i & (n - 1)] = (x[i & (n - 1)] + a[i]) % P;
        for (size_t i = 0; i < b.size(); i++) y[i & (n - 1)] = (y[i & (n - 1)] + b[i]) % P;

        auto rt = roots<P, ROOT>(n, false);
        forward<P>(x, rt), forward<P>(y, rt);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * y[i] % P;

        rt = roots<P, ROOT>(n, true);
        inverse<P>(x, rt);

        unsigned long long inv = power(n, P - 2, P);
        for (auto& v : x) v = v * inv % P;
        return x;
    }

    /// Cyclic convolution of length n modulo m, for values in [0, m)
    vector<long long> cyclic_mod(const vector<long long>& a, const vector<long long>& b, int n, long long m){
        if (m == MOD){
            auto c = cyclic<MOD, 3>(a, b, n);
            return vector<long long>(c.begin(), c.end());
        }

        const unsigned long long INV1 = power(P1, P2 - 2, P2), INV12 = power((long long)P1 * P2 % P3, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (long long)P1 * P2 % m;
        auto c1 = cyclic<P1, 3>(a, b, n), c2 = cyclic<P2, 3>(a, b, n), c3 = cyclic<P3, 11>(a, b, n);

        /// x1 + x2 P1 + x3 P1 P2 is the true coefficient, reduced term by term below 2^62
        vector<long long> res(n);
        for (int i = 0; i < n; i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
        }

        return res;
    }

    /// a * b modulo m, for values in [0, m)
    vector<long long> multiply(const vector<long long>& a, const vector<long long>& b, long long m){
        if (a.empty() || b.empty()) return {};
        int need = a.size() + b.size() - 1, n = 1;
        while (n < need - 1) n <<= 1;

        /// A product one longer than n wraps only its top coefficient a.back() * b.back() onto index 0
        auto res = cyclic_mod(a, b, n, m);
        res.resize(need);
        if (need > n){
            long long top = a.back() * b.back() % m;
            res[0] = (res[0] - top + m) % m, res[n] = top;
        }

        return res;
    }

    bool is_prime(long long x){
        if (x < 2) return false;
        for (long long d = 2; d * d <= x; d++){
            if (x % d == 0) return false;
        }

        return true;
    }

    /// i! and 1 / i! modulo a prime p > n for i = 0..n
    void factorials(int n, long long p, vector<long long>& fact, vector<long long>& inv_fact){
        fact.assign(n + 1, 1), inv_fact.assign(n + 1, 1);
        for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % p;
        inv_fact[n] = power(fact[n], p - 2, p);
        for (int i = n; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % p;
    }
    /// END COPY stirling_ntt
}

/// Fills b_i = B_i / i! for i in [l, r) from i b_i = sum over k < i of b_k / (i - 1 - k)!, acc[i] collects that sum
void bell_online(int l, int r, long long mod, vector<long long>& b, vector<long long>& acc,
                 const vector<long long>& fact, const vector<long long>& inv_fact){
    if (r - l == 1){
        if (l > 0) b[l] = acc[l] * inv_fact[l] % mod * fact[l - 1] % mod;
        return;
    }

    int mid = (l + r) / 2;
    bell_online(l, mid, mod, b, acc, fact, inv_fact);

    /// b_k / (t - 1 - k)! for k in [l, mid) lands on index t - 1 - l of the product, for every t in [mid, r)
    vector<long long> x(b.begin() + l, b.begin() + mid), y(inv_fact.begin(), inv_fact.begin() + (r - l - 1));
    vector<long long> prod = bell_ntt::multiply(x, y, mod);
    for (int t = mid; t < r; t++){
        acc[t] += prod[t - 1 - l];
        if (acc[t] >= mod) acc[t] -= mod;
    }

    bell_online(mid, r, mod, b, acc, fact, inv_fact);
}

vector<long long> bell_numbers_ntt(int n, long long mod = 998244353){
    assert(0 <= n && n <= (1 << 22) && n < mod && mod <= (1 << 30) && bell_ntt::is_prime(mod));

    vector<long long> fact, inv_fact, b(n + 1, 0), acc(n + 1, 0);
    bell_ntt::factorials(n, mod, fact, inv_fact);
    b[0] = 1;
    bell_online(0, n + 1, mod, b, acc, fact, inv_fact);

    vector<long long> res(n + 1);
    for (int i = 0; i <= n; i++) res[i] = b[i] * fact[i] % mod;
    return res;
}

long long bell_number(long long n, long long p){
    assert(0 <= n && n <= (long long)1e18 && 2 <= p && p <= (1 << 20) && bell_ntt::is_prime(p));

    if (n < p) return bell_numbers_ntt(n, p)[n];
    vector<long long> base = bell_numbers_ntt(p - 1, p);

    /// Touchard: B_(i + p) = B_i + B_(i + 1) mod p, so x^n mod x^p - x - 1 = sum c_i x^i gives B_n = sum c_i B_i
    vector<long long> pw(p, 0);
    pw[0] = 1;
    for (int bit = __lg(n); bit >= 0; bit--){
        vector<long long> sq = bell_ntt::multiply(pw, pw, p);

        /// x^i = x^(i - p + 1) + x^(i - p), from the top so every folded term lands on a degree still to be folded
        for (int i = (int)sq.size() - 1; i >= p; i--){
            sq[i - p + 1] = (sq[i - p + 1] + sq[i]) % p;
            sq[i - p] = (sq[i - p] + sq[i]) % p;
        }
        sq.resize(p);
        pw = sq;

        if (n >> bit & 1){
            long long top = pw[p - 1];
            for (int i = p - 1; i > 0; i--) pw[i] = pw[i - 1];
            pw[0] = top;
            pw[1] = (pw[1] + top) % p;
        }
    }

    long long res = 0;
    for (int i = 0; i < p; i++) res = (res + pw[i] * base[i]) % p;
    return res;
}

vector<long long> bell_numbers(int n, long long m){
    assert(n >= 0 && 1 <= m && m <= (1LL << 62));

    /// Bell triangle: row i has i + 1 entries, starts with the last entry of row i - 1, and B_i is its head
    vector<long long> row(n + 1), next(n + 1), res(n + 1);
    row[0] = res[0] = 1 % m;
    for (int i = 1; i <= n; i++){
        next[0] = row[i - 1];
        for (int j = 1; j <= i; j++){
            long long y = next[j - 1] + row[j - 1];
            next[j] = y >= m ? y - m : y;
        }
        swap(row, next);
        res[i] = row[0];
    }

    return res;
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

const int MOD = 1000000007;

int main(){
    int m, n, i, j;

    read(m, n);
    vector<string> ar(n);
    read(ar);

    map<string, int> cnt;
    for (j = 0; j < m; j++){
        string col;
        for (i = 0; i < n; i++) col += ar[i][j];
        cnt[col]++;
    }

    vector<long long> bell = bell_numbers(m, MOD);
    long long res = 1;
    for (auto it: cnt) res = res * bell[it.second] % MOD;

    write(res);
    return 0;
}
