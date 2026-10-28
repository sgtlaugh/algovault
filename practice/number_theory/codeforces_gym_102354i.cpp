/***
 *
 * Chinese Remainder Theorem
 * Solves the system x % mods[i] = rems[i] for arbitrary positive moduli, co-prime or not
 *
 * Complexity, for k congruences and M the largest modulus:
 *   - CRT: O(k log M)
 *   - Garner g(mods, mod): O(k^2 + k log M)
 *   - g(rems): O(k^2)
 *
 * CRT(rems, mods): the unique solution x in [0, lcm(mods)), -1 if the system is inconsistent, 0 if it is empty
 * Garner g(mods, mod): prepares g(rems) for pairwise co-prime mods, asserted, computing the inverses once
 * g(rems): x % mod for the solution x in [0, product of the mods), allocation free, like one call per NTT coefficient
 *
 * Notes:
 *   - CRT: each mod positive, lcm(mods) < 2^63 asserted, rems any long long
 *   - Garner: each mod and the target mod in [1, 2^63), rems any long long, the product of the mods may exceed 64 bits
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY extended_gcd from code_library/number_theory/extended_gcd.cpp
array<long long, 3> extended_gcd(long long a, long long b){
    assert(a >= 0 && b >= 0);

    long long x = 1, y = 0, x1 = 0, y1 = 1;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }

    return {a, x, y};
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    auto [g, x, y] = extended_gcd(a, m);
    if (g != 1) return -1;
    return x < 0 ? x + m : x;
}
/// END COPY extended_gcd

/// x % m + m would overflow for m above 2^62
long long non_negative_mod(long long x, long long m){
    x %= m;
    return x < 0 ? x + m : x;
}

long long CRT(const vector<long long>& rems, const vector<long long>& mods){
    long long res = 0, lcm = 1;

    for (size_t i = 0; i < mods.size(); i++){
        long long r = non_negative_mod(rems[i], mods[i]);
        long long g = gcd(lcm, mods[i]), step = mods[i] / g;
        if ((r - res) % g) return -1;
        assert((__int128)lcm * step <= LLONG_MAX);

        /// res + k * lcm hits r modulo mods[i] exactly when k * (lcm / g) = (r - res) / g modulo step
        long long k = (__int128)non_negative_mod((r - res) / g, step) * mod_inverse(lcm / g, step) % step;
        res += k * lcm;
        lcm *= step;
    }

    return res;
}

/// Garner's digits: x = digit[0] + digit[1] * mods[0] + digit[2] * mods[0] * mods[1] + ..., m[k] is the target mod
struct Garner{
    int k;
    vector<long long> m, inv, coef, digit;

    /// coef[i * (k + 1) + j] = mods[0] * ... * mods[j - 1] % m[i], inv[i] = (mods[0] * ... * mods[i - 1])^-1 % m[i]
    Garner(const vector<long long>& mods, long long mod) : k(mods.size()), m(mods), inv(k), coef((k + 1) * (k + 1)), digit(k){
        m.push_back(mod);

        for (int i = 0; i <= k; i++){
            long long* c = &coef[i * (k + 1)];
            c[0] = 1;
            for (int j = 1; j <= i; j++) c[j] = (unsigned __int128)c[j - 1] * m[j - 1] % m[i];
            if (i < k){
                inv[i] = mod_inverse(c[i], m[i]);
                assert(inv[i] != -1);
            }
        }
    }

    long long operator()(const vector<long long>& rems){
        for (int i = 0; i < k; i++){
            long long diff = non_negative_mod(rems[i], m[i]) - prefix(i);
            if (diff < 0) diff += m[i];
            digit[i] = (unsigned __int128)diff * inv[i] % m[i];
        }

        return prefix(k);
    }

    /// The number written by the first i digits, % m[i]
    long long prefix(int i) const{
        const long long* c = &coef[i * (k + 1)];
        unsigned __int128 res = 0;
        for (int j = 0; j < i; j++){
            res += (unsigned __int128)digit[j] * c[j];
            if (j % 3 == 2) res %= m[i];  /// 3 products below 2^126 plus a residue fit 128 bits, 7-36% faster
        }

        return res % m[i];
    }
};

/***
 *
 * Continued Fractions
 * Continued fraction expansion, convergents, best rational approximation and fraction binary search
 *
 * Complexity: O(log q) for the expansion, convergents and best approximation, O(log N) predicate calls for the search
 *
 * best_approximation(p, q, N): the fraction closest to p / q among all fractions with denominator <= N
 *     Ties go to the smaller denominator, then to the smaller value
 *     Requires p > LLONG_MIN, q >= 1, N >= 1 and __int128
 *
 * continued_fraction(p, q): [a0; a1, ..., ak] with p / q = a0 + 1 / (a1 + 1 / (... + 1 / ak))
 *     a0 = floor(p / q) can be negative, a1..ak >= 1 and ak >= 2 when k >= 1, so the expansion is unique
 *     p / q need not be reduced, requires q >= 1, any p works
 *
 * convergents(a): the fractions [a0; a1, ..., ai] for every i, reduced, the last one is p / q reduced
 *     Requires a non-empty expansion, never overflows on continued_fraction(p, q) for p > LLONG_MIN
 *
 * fraction_search(f, N): the smallest fraction p / q with 0 <= p <= N and 1 <= q <= N for which f(p, q) is true
 *     f must be monotone over the non-negative rationals: false up to some point, true from there on
 *     Returns {1, 0} if f is false for every such fraction
 *     f is only called with reduced fractions inside the bounds, at most 4 log2(N) + 10 times, N <= 1e18
 *
 * Notes:
 *   - Fractions are {numerator, denominator} pairs, the denominator is always positive
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> continued_fraction(long long p, long long q){
    assert(q > 0);
    vector<long long> a;

    while (q){
        long long d = p / q, r = p % q;
        if (r < 0) d--, r += q;
        a.push_back(d);
        p = q, q = r;
    }

    return a;
}

pair<long long, long long> best_approximation(long long p, long long q, long long n){
    assert(q > 0 && n > 0);
    auto a = continued_fraction(p, q);
    long long pp = 1, pq = 0, cp = a[0], cq = 1;

    for (size_t i = 1; i < a.size(); i++){
        if (a[i] <= (n - pq) / cq){
            tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
            continue;
        }

        /// The convergent and the largest semiconvergent are the neighbours of p / q among denominators <= n
        long long t = (n - pq) / cq, sp = pp + t * cp, sq = pq + t * cq;
        __int128 dc = (__int128)p * cq - (__int128)cp * q, ds = (__int128)p * sq - (__int128)sp * q;
        if (dc < 0) dc = -dc;
        if (ds < 0) ds = -ds;

        __int128 lhs = dc * sq, rhs = ds * cq;
        if (lhs != rhs) return lhs < rhs ? make_pair(cp, cq) : make_pair(sp, sq);
        if (cq != sq) return cq < sq ? make_pair(cp, cq) : make_pair(sp, sq);
        return (__int128)cp * sq < (__int128)sp * cq ? make_pair(cp, cq) : make_pair(sp, sq);
    }

    return {cp, cq};
}

vector<pair<long long, long long>> convergents(const vector<long long>& a){
    assert(!a.empty());
    vector<pair<long long, long long>> res;
    long long pp = 1, pq = 0, cp = a[0], cq = 1;
    res.push_back({cp, cq});

    for (size_t i = 1; i < a.size(); i++){
        tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
        res.push_back({cp, cq});
    }

    return res;
}

template<typename F>
pair<long long, long long> fraction_search(F f, long long n){
    assert(n > 0);
    if (f(0LL, 1LL)) return {0, 1};

    /// Invariant: f(lp / lq) is false, f(rp / rq) is true, with 1 / 0 standing for infinity
    long long lp = 0, lq = 1, rp = 1, rq = 0;
    while (lp + rp <= n && lq + rq <= n){
        bool side = f(lp + rp, lq + rq);
        long long& bp = side ? rp : lp;
        long long& bq = side ? rq : lq;
        long long dp = side ? lp : rp, dq = side ? lq : rq;

        auto valid = [&](long long k){
            if (dp && k > (n - bp) / dp) return false;
            if (dq && k > (n - bq) / dq) return false;
            return f(bp + k * dp, bq + k * dq) == side;
        };

        long long k = 1, step = 1;
        while (valid(k + step)) k += step, step *= 2;
        for (step /= 2; step; step /= 2){
            if (valid(k + step)) k += step;
        }

        bp += k * dp, bq += k * dq;
    }

    return {rp, rq};
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

const long long P1 = 2000000011, P2 = 3000000019, LIM = 1000000000;

pair<long long, long long> solve(long long y){
    long long m = P1 * P2;

    for (auto c: convergents(continued_fraction(y, m))){
        if (c.second > LIM) break;
        __int128 p = (__int128)y * c.second - (__int128)c.first * m;
        if (p >= 1 && p <= LIM) return {(long long)p, c.second};
    }

    return {-1, -1};
}

int main(){
    int t;
    long long y1, y2;

    read(t);
    while (t--){
        write('?', P1);
        flush();
        read(y1);

        write('?', P2);
        flush();
        read(y2);

        pair<long long, long long> res = solve(CRT({y1, y2}, {P1, P2}));
        write('!', res.first, res.second);
        flush();
    }

    return 0;
}
