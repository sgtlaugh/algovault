/***
 *
 * Polynomial and Formal Power Series
 * Polynomial arithmetic and power series operations modulo an NTT prime
 *
 * Complexity:
 *   - O(n log n) for *, /, %, divmod, inverse, log, exp, sqrt, pow and taylor_shift
 *   - O(n log^2 n) for evaluate at many points and interpolate
 *   - O(n) for +, -, derivative, integral and evaluate at one point
 *
 * Polynomial<MOD> p(v): v[0] + v[1] x + v[2] x^2 + ... from a vector or braced list, default MOD = 998244353
 * p + q, p - q, p * q: the sum, difference and product
 * p / q, p % q: the quotient and remainder, q must be nonzero
 * p.coef(i): the coefficient of x^i, 0 when i is out of range
 * p.derivative(): the derivative, max(size() - 1, 0) coefficients
 * p.divmod(q): the pair {p / q, p % q}, q must be nonzero
 * p.evaluate(x): the value p(x)
 * p.evaluate(xs): the value p(x) for every point x of xs
 * p.exp(n): e^p mod x^n, requires p.coef(0) = 0
 * p.integral(): the integral with constant term 0, size() + 1 coefficients
 * Polynomial<MOD>::interpolate(xs, ys): the degree < n polynomial through (xs[i], ys[i]), the n xs distinct mod MOD
 * p.inverse(n): 1 / p mod x^n, requires p.coef(0) != 0
 * p.log(n): ln p mod x^n, requires p.coef(0) = 1
 * p.pow(k, n): p^k mod x^n for any k >= 0, p may start with zeros, and 0^0 = 1
 * p.size(): the number of coefficients, trailing zeros included
 * p.sqrt(n): a square root of p mod x^n, lowest nonzero coefficient the smaller modular root, nullopt when none exists
 * p.taylor_shift(c): p(x + c), size() coefficients
 * p.trimmed(): p without trailing zeros
 * p.a: the coefficients in [0, MOD), lowest degree first
 *
 * Notes:
 *   - MOD is a prime below 2^30 of the form c * 2^v + 1, the default 998244353 = 119 * 2^23 + 1
 *   - Every series length n, polynomial length and number of points is at most 2^(v - 1), 2^22 for 998244353
 *   - Results of +, -, *, /, %, divmod and interpolate have trailing zeros trimmed, and the zero polynomial is empty
 *   - inverse, log, exp, sqrt and pow return exactly n coefficients and read only the first n of p
 *   - Coefficients, points and c may be any integers, they are reduced into [0, MOD)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<int MOD = 998244353>
struct Polynomial{
    static_assert(MOD > 2 && MOD < (1 << 30), "MOD must be an odd NTT prime below 2^30");

    vector<int> a;

    Polynomial() {}

    Polynomial(initializer_list<long long> v){
        for (auto x : v) a.push_back(norm(x));
    }

    template<typename T>
    Polynomial(const vector<T>& v){
        a.reserve(v.size());
        for (auto x : v) a.push_back(norm(x));
    }

    Polynomial operator + (const Polynomial& b) const{
        vector<int> res(max(size(), b.size()));
        for (int i = 0; i < (int)res.size(); i++) res[i] = add(coef(i), b.coef(i));
        return Polynomial(res).trimmed();
    }

    Polynomial operator - (const Polynomial& b) const{
        vector<int> res(max(size(), b.size()));
        for (int i = 0; i < (int)res.size(); i++) res[i] = sub(coef(i), b.coef(i));
        return Polynomial(res).trimmed();
    }

    Polynomial operator * (const Polynomial& b) const{
        return Polynomial(multiply(a, b.a)).trimmed();
    }

    Polynomial operator / (const Polynomial& b) const{
        return divmod(b).first;
    }

    Polynomial operator % (const Polynomial& b) const{
        return divmod(b).second;
    }

    int coef(int i) const{
        return (i >= 0 && i < size()) ? a[i] : 0;
    }

    Polynomial derivative() const{
        vector<int> res(max(size() - 1, 0));
        for (int i = 1; i < size(); i++) res[i - 1] = mul(a[i], i);
        return res;
    }

    pair<Polynomial, Polynomial> divmod(const Polynomial& b) const{
        Polynomial p = trimmed(), q = b.trimmed();

        assert(q.size() > 0);
        if (p.size() < q.size()) return {Polynomial(), p};

        int d = p.size() - q.size() + 1;
        vector<int> rp(p.a.rbegin(), p.a.rbegin() + d), rq(q.a.rbegin(), q.a.rbegin() + min(d, q.size()));
        vector<int> quot = multiply(rp, Polynomial(rq).inverse(d).a);
        quot.resize(d);
        reverse(quot.begin(), quot.end());

        Polynomial res(quot);
        return {res, p - res * q};
    }

    int evaluate(long long x) const{
        int v = norm(x), res = 0;
        for (int i = size() - 1; i >= 0; i--) res = add(mul(res, v), a[i]);
        return res;
    }

    vector<int> evaluate(const vector<int>& xs) const{
        int n = xs.size();
        vector<int> pts(n), res(n);

        for (int i = 0; i < n; i++) pts[i] = norm(xs[i]);
        if (n == 0) return res;

        vector<Polynomial> tree(4 * n);
        build(tree, 1, 0, n - 1, pts);
        descend(tree, 1, 0, n - 1, *this, pts, res);
        return res;
    }

    Polynomial exp(int n) const{
        assert(n >= 0);
        if (n == 0) return {};
        assert(coef(0) == 0);

        /// Newton step g = g * (1 + P - ln g) doubles the number of correct coefficients
        Polynomial g = {1};
        for (int m = 1; m < n; m *= 2){
            vector<int> h = prefix(2 * m), lg = g.log(2 * m).a;
            for (int i = 0; i < 2 * m; i++) h[i] = sub(h[i], lg[i]);
            h[0] = add(h[0], 1);
            g.a = multiply(g.a, h);
            g.a.resize(2 * m);
        }

        g.a.resize(n);
        return g;
    }

    Polynomial integral() const{
        int n = size();
        vector<int> inv(n + 1, 1), res(n + 1);
        for (int i = 2; i <= n; i++) inv[i] = mul(MOD - MOD / i, inv[MOD % i]);
        for (int i = 0; i < n; i++) res[i + 1] = mul(a[i], inv[i + 1]);
        return res;
    }

    static Polynomial interpolate(const vector<int>& xs, const vector<int>& ys){
        assert(xs.size() == ys.size());
        int n = xs.size();
        if (n == 0) return {};

        vector<int> pts(n), w(n);
        for (int i = 0; i < n; i++) pts[i] = norm(xs[i]);
        vector<Polynomial> tree(4 * n);
        build(tree, 1, 0, n - 1, pts);

        /// w[i] = y[i] / M'(x[i]) for M = prod (x - x[j]), M'(x[i]) is zero exactly when a point repeats
        descend(tree, 1, 0, n - 1, tree[1].derivative(), pts, w);
        for (int i = 0; i < n; i++){
            assert(w[i] != 0);
            w[i] = mul(norm(ys[i]), power(w[i], MOD - 2));
        }

        return combine(tree, 1, 0, n - 1, w).trimmed();
    }

    Polynomial inverse(int n) const{
        assert(n >= 0);
        if (n == 0) return {};
        assert(coef(0) != 0);

        vector<int> b = {power(a[0], MOD - 2)};
        for (int m = 1; m < n; m *= 2){
            vector<int> f = prefix(2 * m), g = b, roots = make_roots(2 * m);
            g.resize(2 * m);
            ntt(f, roots);
            ntt(g, roots);
            for (int i = 0; i < 2 * m; i++) f[i] = mul(f[i], g[i]);
            inverse_ntt(f, roots);

            /// Cyclic products of length 2m: the wrapped terms land below m, [m, 2m) is exact, and P * b = 1 below m
            fill(f.begin(), f.begin() + m, 0);
            ntt(f, roots);
            for (int i = 0; i < 2 * m; i++) f[i] = mul(f[i], g[i]);
            inverse_ntt(f, roots);

            b.resize(2 * m);
            for (int i = m; i < 2 * m; i++) b[i] = sub(0, f[i]);
        }

        b.resize(n);
        return b;
    }

    Polynomial log(int n) const{
        assert(n >= 0);
        if (n == 0) return {};
        assert(coef(0) == 1);

        vector<int> q = multiply(Polynomial(prefix(n)).derivative().a, inverse(n).a);
        q.resize(n - 1);
        return Polynomial(q).integral();
    }

    Polynomial pow(long long k, int n) const{
        assert(k >= 0 && n >= 0);
        vector<int> res(n);
        if (n == 0) return res;
        if (k == 0){
            res[0] = 1;
            return res;
        }

        int i = 0, lim = min(n, size());
        while (i < lim && a[i] == 0) i++;
        if (i == lim || (i > 0 && k >= (n + i - 1) / i)) return res;

        /// P = c x^i T with T[0] = 1, T^k = exp(k ln T) where k only matters modulo MOD
        int shift = i * k, len = n - shift, c = a[i], c_inv = power(c, MOD - 2);
        vector<int> t(len);
        for (int j = 0; j < len && i + j < size(); j++) t[j] = mul(a[i + j], c_inv);
        Polynomial lg = Polynomial(t).log(len);
        for (auto& x : lg.a) x = mul(x, k % MOD);

        Polynomial e = lg.exp(len);
        int ck = power(c, k);
        for (int j = 0; j < len; j++) res[shift + j] = mul(e.a[j], ck);
        return res;
    }

    int size() const{
        return a.size();
    }

    Polynomial taylor_shift(long long c) const{
        int n = size(), v = norm(c);

        if (n == 0) return {};

        vector<int> fact(n, 1), inv_fact(n);
        for (int i = 1; i < n; i++) fact[i] = mul(fact[i - 1], i);
        inv_fact[n - 1] = power(fact[n - 1], MOD - 2);
        for (int i = n - 1; i > 0; i--) inv_fact[i - 1] = mul(inv_fact[i], i);

        /// res[j] j! = sum over i >= j of (a[i] i!) (c^(i - j) / (i - j)!), a convolution after reversing a[i] i!
        vector<int> x(n), y(n);
        int pw = 1;
        for (int i = 0; i < n; i++){
            x[n - 1 - i] = mul(a[i], fact[i]);
            y[i] = mul(pw, inv_fact[i]);
            pw = mul(pw, v);
        }

        vector<int> z = multiply(x, y), res(n);
        for (int j = 0; j < n; j++) res[j] = mul(z[n - 1 - j], inv_fact[j]);
        return res;
    }

    Polynomial trimmed() const{
        Polynomial res = *this;
        while (!res.a.empty() && res.a.back() == 0) res.a.pop_back();
        return res;
    }

    static int norm(long long x){
        x %= MOD;
        return x < 0 ? x + MOD : x;
    }

    static int add(int x, int y){
        return x + y >= MOD ? x + y - MOD : x + y;
    }

    static int sub(int x, int y){
        return x - y < 0 ? x - y + MOD : x - y;
    }

    static vector<int> multiply(const vector<int>& x, const vector<int>& y){
        if (x.empty() || y.empty()) return {};

        int n = x.size(), m = y.size();
        vector<int> res(n + m - 1);

        /// Without this schoolbook cutoff interpolate at 2^17 points measured 12% slower
        if (min(n, m) <= 32){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++) res[i + j] = add(res[i + j], mul(x[i], y[j]));
            }

            return res;
        }

        int len = 1;
        while (len < n + m - 1) len *= 2;
        vector<int> f = x, g = y, roots = make_roots(len);
        f.resize(len);
        g.resize(len);
        ntt(f, roots);
        ntt(g, roots);
        for (int i = 0; i < len; i++) f[i] = mul(f[i], g[i]);
        inverse_ntt(f, roots);

        f.resize(n + m - 1);
        return f;
    }

    static int mul(int x, int y){
        return (long long)x * y % MOD;
    }

    /// tree[v] = prod (x - pts[i]) over the node's range
    static void build(vector<Polynomial>& tree, int v, int l, int r, const vector<int>& pts){
        if (l == r){
            tree[v].a = {sub(0, pts[l]), 1};
            return;
        }

        int mid = (l + r) / 2;
        build(tree, 2 * v, l, mid, pts);
        build(tree, 2 * v + 1, mid + 1, r, pts);
        tree[v].a = multiply(tree[2 * v].a, tree[2 * v + 1].a);
    }

    /// p mod tree[v] agrees with p at every point of the node's range
    static void descend(const vector<Polynomial>& tree, int v, int l, int r, Polynomial p, const vector<int>& pts, vector<int>& res){
        p = p % tree[v];

        /// Horner below 32 points: descending to single points measured 30% slower at 2^17 points
        if (r - l < 32){
            for (int i = l; i <= r; i++) res[i] = p.evaluate(pts[i]);
            return;
        }

        int mid = (l + r) / 2;
        descend(tree, 2 * v, l, mid, p, pts, res);
        descend(tree, 2 * v + 1, mid + 1, r, p, pts, res);
    }

    vector<int> prefix(int k) const{
        vector<int> res(k);
        for (int i = 0; i < min(k, size()); i++) res[i] = a[i];
        return res;
    }

    static int power(long long b, long long e){
        long long res = 1;

        b %= MOD;
        while (e > 0){
            if (e & 1) res = res * b % MOD;
            b = b * b % MOD;
            e >>= 1;
        }

        return res;
    }

    /// Sum of w[i] tree[v] / (x - pts[i]) over the node's range
    static Polynomial combine(const vector<Polynomial>& tree, int v, int l, int r, const vector<int>& w){
        if (l == r) return vector<int>{w[l]};

        int mid = (l + r) / 2;
        return combine(tree, 2 * v, l, mid, w) * tree[2 * v + 1] + combine(tree, 2 * v + 1, mid + 1, r, w) * tree[2 * v];
    }

    /// roots[k + j] = w^j for the principal 2k-th root of unity w, for every power of two k < len
    static vector<int> make_roots(int len){
        assert((MOD - 1) % len == 0);
        int g = primitive_root();

        vector<int> roots(max(len, 2));
        roots[1] = 1;
        for (int k = 1; 2 * k < len; k *= 2){
            int w = power(g, (MOD - 1) / (4 * k));
            for (int j = k; j < 2 * k; j++){
                roots[2 * j] = roots[j];
                roots[2 * j + 1] = mul(roots[j], w);
            }
        }

        return roots;
    }

    static void ntt(vector<int>& f, const vector<int>& roots){
        int n = f.size();

        for (int i = 1, j = 0; i < n; i++){
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(f[i], f[j]);
        }

        for (int k = 1; k < n; k *= 2){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    int z = mul(f[i + j + k], roots[j + k]);
                    f[i + j + k] = sub(f[i + j], z);
                    f[i + j] = add(f[i + j], z);
                }
            }
        }
    }

    static void inverse_ntt(vector<int>& f, const vector<int>& roots){
        ntt(f, roots);
        reverse(f.begin() + 1, f.end());

        int inv = power(f.size(), MOD - 2);
        for (auto& x : f) x = mul(x, inv);
    }

    /// Tonelli-Shanks, -1 when x is not a quadratic residue, otherwise the smaller of the two roots
    static int mod_sqrt(int x){
        if (x == 0) return 0;
        if (power(x, (MOD - 1) / 2) != 1) return -1;

        int s = 0, q = MOD - 1, z = 2;
        while (q % 2 == 0){
            q /= 2;
            s++;
        }
        while (power(z, (MOD - 1) / 2) != MOD - 1) z++;

        int c = power(z, q), t = power(x, q), r = power(x, (q + 1) / 2);
        while (t != 1){
            int i = 0;
            for (int u = t; u != 1; u = mul(u, u)) i++;

            int b = c;
            for (int j = 0; j < s - i - 1; j++) b = mul(b, b);
            s = i;
            c = mul(b, b);
            t = mul(t, c);
            r = mul(r, b);
        }

        return min(r, MOD - r);
    }

    static int primitive_root(){
        int factors[32] = {}, cnt = 0, x = MOD - 1;

        for (int d = 2; d * d <= x; d++){
            if (x % d) continue;
            factors[cnt++] = d;
            while (x % d == 0) x /= d;
        }
        if (x > 1) factors[cnt++] = x;

        for (int g = 2; ; g++){
            bool ok = true;
            for (int i = 0; i < cnt && ok; i++) ok = power(g, (MOD - 1) / factors[i]) != 1;
            if (ok) return g;
        }
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

template <int MOD>
bool squarefree(const vector<int>& c){
    Polynomial<MOD> p(c), q = p.derivative();
    while (q.size()){
        Polynomial<MOD> r = p % q;
        p = q, q = r;
    }

    return p.size() == 1;
}

bool distinct(const vector<int>& c){
    return squarefree<998244353>(c) || squarefree<469762049>(c) || squarefree<167772161>(c)
        || squarefree<754974721>(c) || squarefree<1004535809>(c) || squarefree<985661441>(c);
}

int main(){
    int t, n, i;

    while (read(t)){
        while (t--){
            read(n);
            vector<int> c(n + 1);
            for (i = n; i >= 0; i--) read(c[i]);

            if (n > 0 && c[n] == 0){
                write("No!");
                continue;
            }

            write(distinct(c) ? "Yes!" : "No!");
        }
    }

    return 0;
}
