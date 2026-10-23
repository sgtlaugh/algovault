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
 * Polynomial<MOD> p(v): the polynomial v[0] + v[1] x + v[2] x^2 + ... from a vector or braced list of integers, default MOD = 998244353
 * p.a: the coefficients in [0, MOD), lowest degree first
 * p.size(): the number of coefficients, trailing zeros included
 * p.coef(i): the coefficient of x^i, 0 when i is out of range
 * p.trimmed(): p without trailing zeros
 * p + q, p - q, p * q: the sum, difference and product
 * p / q, p % q: the quotient and remainder, q must be nonzero
 * p.divmod(q): the pair {p / q, p % q}, q must be nonzero
 * p.derivative(): the derivative, max(size() - 1, 0) coefficients
 * p.integral(): the integral with constant term 0, size() + 1 coefficients
 * p.evaluate(x): the value p(x)
 * p.evaluate(xs): the value p(x) for every point x of xs
 * Polynomial<MOD>::interpolate(xs, ys): the polynomial of degree < n through the n points (xs[i], ys[i]), the xs distinct modulo MOD
 * p.inverse(n): 1 / p mod x^n, requires p.coef(0) != 0
 * p.log(n): ln p mod x^n, requires p.coef(0) = 1
 * p.exp(n): e^p mod x^n, requires p.coef(0) = 0
 * p.sqrt(n): a square root of p mod x^n whose lowest nonzero coefficient is the smaller modular root, nullopt when none exists
 * p.pow(k, n): p^k mod x^n for any k >= 0, p may start with zeros, and 0^0 = 1
 * p.taylor_shift(c): p(x + c), size() coefficients
 *
 * MOD is a prime below 2^30 of the form c * 2^v + 1, the default 998244353 = 119 * 2^23 + 1
 * Every size involved, the series length n, the polynomial lengths and the number of points, must be at most 2^(v - 1), 2^22 for 998244353
 * Results of +, -, *, /, %, divmod and interpolate have trailing zeros trimmed, and the zero polynomial is empty
 * inverse, log, exp, sqrt and pow return exactly n coefficients and read only the first n of p
 * Coefficients, points and c may be any integers, they are reduced into [0, MOD)
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

    int size() const{
        return a.size();
    }

    int coef(int i) const{
        return (i >= 0 && i < size()) ? a[i] : 0;
    }

    Polynomial trimmed() const{
        Polynomial res = *this;
        while (!res.a.empty() && res.a.back() == 0) res.a.pop_back();
        return res;
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

    Polynomial derivative() const{
        vector<int> res(max(size() - 1, 0));
        for (int i = 1; i < size(); i++) res[i - 1] = mul(a[i], i);
        return res;
    }

    Polynomial integral() const{
        int n = size();
        vector<int> inv(n + 1, 1), res(n + 1);
        for (int i = 2; i <= n; i++) inv[i] = mul(MOD - MOD / i, inv[MOD % i]);
        for (int i = 0; i < n; i++) res[i + 1] = mul(a[i], inv[i + 1]);
        return res;
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

    optional<Polynomial> sqrt(int n) const{
        assert(n >= 0);
        vector<int> res(n);
        int i = 0, lim = min(n, size());
        while (i < lim && a[i] == 0) i++;
        if (i == lim) return Polynomial(res);
        if (i % 2) return nullopt;

        int s = mod_sqrt(a[i]);
        if (s < 0) return nullopt;

        /// P mod x^n = x^i Y, the root is x^(i / 2) R with R^2 = Y, Newton step R = (R + Y / R) / 2
        int len = n - i / 2, inv2 = (MOD + 1) / 2;
        Polynomial y(vector<int>(a.begin() + i, a.begin() + lim)), r = {s};
        for (int m = 1; m < len; m *= 2){
            vector<int> q = multiply(y.prefix(2 * m), r.inverse(2 * m).a);
            r.a.resize(2 * m);
            for (int j = 0; j < 2 * m; j++) r.a[j] = mul(add(r.a[j], q[j]), inv2);
        }

        for (int j = 0; j < len; j++) res[i / 2 + j] = r.a[j];
        return Polynomial(res);
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

    static int mul(int x, int y){
        return (long long)x * y % MOD;
    }

    static constexpr int power(long long b, long long e){
        long long res = 1;
        b %= MOD;
        while (e > 0){
            if (e & 1) res = res * b % MOD;
            b = b * b % MOD;
            e >>= 1;
        }
        return res;
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

    /// roots[k + j] = w^j for the principal 2k-th root of unity w, for every power of two k < len
    static vector<int> make_roots(int len){
        assert((MOD - 1) % len == 0);
        constexpr int g = primitive_root();

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

    static constexpr int primitive_root(){
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

    vector<int> prefix(int k) const{
        vector<int> res(k);
        for (int i = 0; i < min(k, size()); i++) res[i] = a[i];
        return res;
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

    /// Sum of w[i] tree[v] / (x - pts[i]) over the node's range
    static Polynomial combine(const vector<Polynomial>& tree, int v, int l, int r, const vector<int>& w){
        if (l == r) return vector<int>{w[l]};

        int mid = (l + r) / 2;
        return combine(tree, 2 * v, l, mid, w) * tree[2 * v + 1] + combine(tree, 2 * v + 1, mid + 1, r, w) * tree[2 * v];
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

/// Because is_integral is not always true for __int128
#ifdef __SIZEOF_INT128__
    namespace std{
        template<> struct is_integral<__int128_t>: true_type {};
    }
#endif

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

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
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

    /// enable_if in the return type, since a second default template argument would redeclare the integral overload
    template <typename T>
    typename enable_if<is_floating_point<T>::value, bool>::type read_one(T& x){
        int c = ' ', len = 0;
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;

        /// Parsing at the target's own width, since strtold then a cast to double rounds twice and can miss by one ulp
        if (is_same<T, float>::value) x = strtof(tmpbuf, nullptr);
        else if (is_same<T, double>::value) x = strtod(tmpbuf, nullptr);
        else x = strtold(tmpbuf, nullptr);
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

    bool read_one(char& ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T>& v){
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

    int read(){
        return 0;
    }

    template <typename T, typename ...Args>
    int read(T &x, Args& ...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(const char& c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char* s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    void write_one(const string& s){
        for (auto &&c: s) write_char(c);
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
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

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

    void write_one(const Fixed& f){
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
    void write_one(const vector<T>& v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    void write() {}

    template <typename T> struct is_vector : false_type {};
    template <typename T, typename A> struct is_vector<vector<T, A>> : true_type {};

    template <typename T, typename ...Args>
    void write(const T& x, const Args& ...args){
        write_one(x);
        write_char(sizeof...(args) && !is_vector<T>::value ? ' ' : '\n');
        write(args...);
    }
}

int main(){
    using namespace fio;
    int n, m;

    read(n, m);
    vector<int> f(n), g(m);
    read(f, g);

    auto res = Polynomial<>(f).divmod(Polynomial<>(g));
    write(res.first.size(), res.second.size());
    write(res.first.a);
    write(res.second.a);
    return 0;
}
