/***
 *
 * Polynomial and Formal Power Series
 * Polynomial arithmetic and power series operations modulo an NTT prime, with its own number theoretic transform
 *
 * Complexity:
 *   - O(n log n) for multiply, inverse, divmod, log, exp, sqrt, pow and taylor_shift
 *   - O(n log^2 n) for multipoint evaluate and interpolate
 *   - O(n) for derivative, integral and single point evaluate
 *
 * Polynomial<MOD> holds a[0] + a[1] x + a[2] x^2 + ..., coefficients and points are reduced into [0, MOD)
 * MOD is a prime below 2^30 of the form c * 2^v + 1, the default 998244353 = 119 * 2^23 + 1
 * Every size involved (series length n, polynomial lengths, number of points) must be at most 2^(v - 1),
 * which is 2^22 for the default modulus
 *
 * Polynomial results (+, -, *, /, %, divmod, interpolate) have trailing zeros trimmed, the zero polynomial is empty
 * Power series results (inverse, log, exp, sqrt, pow) have exactly n coefficients and read only the first n of P
 * derivative, integral and taylor_shift are not trimmed, they return max(size() - 1, 0), size() + 1 and size()
 * coefficients
 *
 *   inverse(n)            1 / P mod x^n, requires P[0] != 0
 *   log(n)                ln P mod x^n, requires P[0] = 1
 *   exp(n)                e^P mod x^n, requires P[0] = 0
 *   sqrt(n)               first n coefficients of the power series square root of (P mod x^n) whose lowest
 *                         nonzero coefficient is the smaller modular root, nullopt when no root exists
 *   pow(k, n)             P^k mod x^n for any k >= 0, P may start with zeros, 0^0 = 1
 *   divmod(B)             quotient and remainder of P / B, B must be nonzero
 *   evaluate(xs)          P(x) for every point x
 *   interpolate(xs, ys)   the polynomial of degree < n through n points, the xs must be distinct
 *   taylor_shift(c)       P(x + c)
 *
 * Example:
 *   Polynomial<> p = {0, 1};
 *   p.exp(4).a;      // {1, 1, 1/2, 1/6} = {1, 1, 499122177, 166374059}
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
            ntt(f, roots), ntt(g, roots);
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

    Polynomial operator % (const Polynomial& b) const{
        return divmod(b).second;
    }

    Polynomial operator * (const Polynomial& b) const{
        return Polynomial(multiply(a, b.a)).trimmed();
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

    Polynomial operator / (const Polynomial& b) const{
        return divmod(b).first;
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

    Polynomial taylor_shift(long long c) const{
        int n = size(), v = norm(c);
        if (n == 0) return {};

        vector<int> fact(n, 1), inv_fact(n);
        for (int i = 1; i < n; i++) fact[i] = mul(fact[i - 1], i);
        inv_fact[n - 1] = power(fact[n - 1], MOD - 2);
        for (int i = n - 1; i > 0; i--) inv_fact[i - 1] = mul(inv_fact[i], i);

        /// res[j] j! = sum over i >= j of (a[i] i!) (c^(i - j) / (i - j)!), a convolution after reversing a[i] i!
        vector<int> x(n), y(n);
        for (int i = 0, pw = 1; i < n; i++, pw = mul(pw, v)){
            x[n - 1 - i] = mul(a[i], fact[i]);
            y[i] = mul(pw, inv_fact[i]);
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

private:
    static int add(int x, int y){
        return x + y >= MOD ? x + y - MOD : x + y;
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

    static Polynomial combine(const vector<Polynomial>& tree, int v, int l, int r, const vector<int>& w){
        if (l == r) return vector<int>{w[l]};

        int mid = (l + r) / 2;
        return combine(tree, 2 * v, l, mid, w) * tree[2 * v + 1] + combine(tree, 2 * v + 1, mid + 1, r, w) * tree[2 * v];
    }

    static void descend(const vector<Polynomial>& tree, int v, int l, int r, Polynomial p, const vector<int>& pts, vector<int>& res){
        p = p % tree[v];
        if (r - l < 32){
            for (int i = l; i <= r; i++) res[i] = p.evaluate(pts[i]);
            return;
        }

        int mid = (l + r) / 2;
        descend(tree, 2 * v, l, mid, p, pts, res);
        descend(tree, 2 * v + 1, mid + 1, r, p, pts, res);
    }

    static void inverse_ntt(vector<int>& f, const vector<int>& roots){
        ntt(f, roots);
        reverse(f.begin() + 1, f.end());

        int inv = power(f.size(), MOD - 2);
        for (auto& x : f) x = mul(x, inv);
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

    /// Tonelli-Shanks, -1 when x is not a quadratic residue, otherwise the smaller of the two roots
    static int mod_sqrt(int x){
        if (x == 0) return 0;
        if (power(x, (MOD - 1) / 2) != 1) return -1;

        int s = 0, q = MOD - 1, z = 2;
        while (q % 2 == 0) q /= 2, s++;
        while (power(z, (MOD - 1) / 2) != MOD - 1) z++;

        int c = power(z, q), t = power(x, q), r = power(x, (q + 1) / 2);
        while (t != 1){
            int i = 0;
            for (int u = t; u != 1; u = mul(u, u)) i++;

            int b = c;
            for (int j = 0; j < s - i - 1; j++) b = mul(b, b);
            s = i, c = mul(b, b), t = mul(t, c), r = mul(r, b);
        }

        return min(r, MOD - r);
    }

    static int mul(int x, int y){
        return (long long)x * y % MOD;
    }

    static vector<int> multiply(const vector<int>& x, const vector<int>& y){
        if (x.empty() || y.empty()) return {};

        int n = x.size(), m = y.size();
        vector<int> res(n + m - 1);
        if (min(n, m) <= 32){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++) res[i + j] = add(res[i + j], mul(x[i], y[j]));
            }
            return res;
        }

        int len = 1;
        while (len < n + m - 1) len *= 2;
        vector<int> f = x, g = y, roots = make_roots(len);
        f.resize(len), g.resize(len);
        ntt(f, roots), ntt(g, roots);
        for (int i = 0; i < len; i++) f[i] = mul(f[i], g[i]);
        inverse_ntt(f, roots);

        f.resize(n + m - 1);
        return f;
    }

    static int norm(long long x){
        x %= MOD;
        return x < 0 ? x + MOD : x;
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

    static constexpr int power(long long b, long long e){
        long long res = 1;
        b %= MOD;
        for (; e > 0; e >>= 1, b = b * b % MOD){
            if (e & 1) res = res * b % MOD;
        }
        return res;
    }

    vector<int> prefix(int k) const{
        vector<int> res(k);
        for (int i = 0; i < min(k, size()); i++) res[i] = a[i];
        return res;
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

    static int sub(int x, int y){
        return x - y < 0 ? x - y + MOD : x - y;
    }
};

int main(){
    using P = Polynomial<>;
    const int mod = 998244353, inv2 = 499122177, inv3 = 332748118, inv4 = 748683265, inv6 = 166374059;

    assert((P({-1, mod, 2LL * mod + 5}).a == vector<int>{mod - 1, 0, 5}));
    assert(((P{1, 2, 3} * P{4, 5}).a == vector<int>{4, 13, 22, 15}));
    assert(((P{1, 2, 3} + P{-1, -2, -3}).a == vector<int>{}));
    assert(((P{1, 2} - P{0, 0, 7}).a == vector<int>{1, 2, mod - 7}));

    assert((P{1, -1}.inverse(5).a == vector<int>{1, 1, 1, 1, 1}));
    assert((P{1, 1}.inverse(4).a == vector<int>{1, mod - 1, 1, mod - 1}));
    assert((P{2}.inverse(2).a == vector<int>{inv2, 0}));

    assert(((P{-1, 0, 0, 1} / P{-1, 1}).a == vector<int>{1, 1, 1}));
    assert(((P{-1, 0, 0, 1} % P{-1, 1}).a == vector<int>{}));
    assert(((P{1, 0, 1} % P{-1, 1}).a == vector<int>{2}));
    assert(((P{1, 2} / P{0, 0, 1}).a == vector<int>{}));
    assert(((P{6, 4} / P{2}).a == vector<int>{3, 2}));

    assert((P{1, 2, 3}.derivative().a == vector<int>{2, 6}));
    assert((P{2, 6}.integral().a == vector<int>{0, 2, 3}));
    assert((P{1, 0}.derivative().a == vector<int>{0}));
    assert((P{}.derivative().a == vector<int>{}));
    assert((P{0}.integral().a == vector<int>{0, 0}));
    assert((P{}.integral().a == vector<int>{0}));

    assert((P{1, 1, 1, 1, 1}.log(5).a == vector<int>{0, 1, inv2, inv3, inv4}));
    assert((P{0, 1}.exp(4).a == vector<int>{1, 1, inv2, inv6}));
    assert((P{}.exp(3).a == vector<int>{1, 0, 0}));
    assert((Polynomial<17>{0, 1}.exp(4).a == vector<int>{1, 1, 9, 3}));

    assert((P{1, 1}.sqrt(3)->a == vector<int>{1, inv2, 124780544}));
    assert((P{0, 0, 4}.sqrt(3)->a == vector<int>{0, 2, 0}));
    assert((P{9, 0}.sqrt(2)->a == vector<int>{3, 0}));
    assert((!P{0, 1}.sqrt(3).has_value()));
    assert((!P{3}.sqrt(2).has_value()));
    assert((P{0, 0, 0, 1}.sqrt(2)->a == vector<int>{0, 0}));

    assert((P{0, 1, 1}.pow(3, 6).a == vector<int>{0, 0, 0, 1, 3, 3}));
    assert((P{0, 1}.pow(1000000000000000000LL, 4).a == vector<int>{0, 0, 0, 0}));
    assert((P{0, 1}.pow(0, 3).a == vector<int>{1, 0, 0}));
    assert((P{1, 1}.pow(mod, 3).a == vector<int>{1, 0, 0}));
    assert((P{2, 2}.pow(2, 3).a == vector<int>{4, 8, 4}));

    assert((P{1, 2, 3}.evaluate(vector<int>{0, 1, 2, -1}) == vector<int>{1, 6, 17, 2}));
    assert((P{1, 2, 3}.evaluate(-1) == 2));
    assert((P::interpolate({0, 1, 2}, {1, 6, 17}).a == vector<int>{1, 2, 3}));
    assert((P::interpolate({5}, {-3}).a == vector<int>{mod - 3}));
    assert((P::interpolate({4, 7}, {0, 0}).a == vector<int>{}));
    assert((P::interpolate({4}, {0}).a == vector<int>{}));

    assert((P{0, 0, 1}.taylor_shift(1).a == vector<int>{1, 2, 1}));
    assert((P{1, 2, 3}.taylor_shift(-1).a == vector<int>{2, mod - 4, 3}));
    assert((P{1, 0}.taylor_shift(2).a == vector<int>{1, 0}));

    /// Identities on a series long enough for every transform path: P * P^-1 = 1, exp(log P) = P, sqrt(P^2) = P
    const int n = 5000;
    vector<long long> coefficients(n);
    for (int i = 0; i < n; i++) coefficients[i] = (1LL * i * i * 7919 + 13) % mod;
    coefficients[0] = 1;
    P p(coefficients);

    vector<int> one(n);
    one[0] = 1;
    auto product = (p * p.inverse(n)).a;
    product.resize(n);
    assert(product == one);
    assert(p.log(n).exp(n).a == p.a);
    assert((p * p).sqrt(n)->a == p.a);

    return 0;
}
