#include "../common.h"

#define main library_main
#include "../../code_library/algebra/polynomial.cpp"
#undef main

typedef long long ll;

/// O(n^2) schoolbook versions of every operation, each written from its defining recurrence
template<int MOD>
struct Schoolbook{
    static ll pw(ll b, ll e){
        ll res = 1;
        for (b %= MOD; e > 0; e >>= 1, b = b * b % MOD) if (e & 1) res = res * b % MOD;
        return res;
    }

    static ll inv(ll x){
        return pw(x, MOD - 2);
    }

    static vector<int> trim(vector<int> v){
        while (!v.empty() && v.back() == 0) v.pop_back();
        return v;
    }

    static vector<int> cut(vector<int> v, int n){
        v.resize(n, 0);
        return v;
    }

    static vector<int> multiply(const vector<int>& x, const vector<int>& y){
        if (x.empty() || y.empty()) return {};
        vector<int> res(x.size() + y.size() - 1, 0);
        for (size_t i = 0; i < x.size(); i++){
            for (size_t j = 0; j < y.size(); j++) res[i + j] = (res[i + j] + (ll)x[i] * y[j]) % MOD;
        }
        return res;
    }

    static vector<int> add(const vector<int>& x, const vector<int>& y, ll sign){
        vector<int> res(max(x.size(), y.size()), 0);
        for (size_t i = 0; i < res.size(); i++){
            ll u = i < x.size() ? x[i] : 0, v = i < y.size() ? y[i] : 0;
            res[i] = ((u + sign * v) % MOD + MOD) % MOD;
        }
        return trim(res);
    }

    /// b[0] = 1 / a[0], b[k] = -(a[1] b[k - 1] + ... + a[k] b[0]) / a[0]
    static vector<int> inverse(vector<int> a, int n){
        if (n == 0) return {};
        a = cut(a, n);
        vector<int> b(n, 0);
        ll i0 = inv(a[0]);
        for (int k = 0; k < n; k++){
            ll s = k == 0 ? 1 : 0;
            for (int j = 1; j <= k; j++) s = (s - (ll)a[j] * b[k - j]) % MOD;
            b[k] = (s % MOD + MOD) % MOD * i0 % MOD;
        }
        return b;
    }

    static pair<vector<int>, vector<int>> divmod(vector<int> a, vector<int> b){
        a = trim(a), b = trim(b);
        if (a.size() < b.size()) return {{}, a};

        vector<int> q(a.size() - b.size() + 1, 0);
        ll lead = inv(b.back());
        for (int i = q.size() - 1; i >= 0; i--){
            ll c = a[i + b.size() - 1] * lead % MOD;
            q[i] = c;
            for (size_t j = 0; j < b.size(); j++) a[i + j] = ((a[i + j] - c * b[j]) % MOD + MOD) % MOD;
        }
        return {trim(q), trim(a)};
    }

    static vector<int> derivative(const vector<int>& a){
        vector<int> res;
        for (size_t i = 1; i < a.size(); i++) res.push_back((ll)a[i] * i % MOD);
        return res;
    }

    static vector<int> integral(const vector<int>& a){
        vector<int> res(a.size() + 1, 0);
        for (size_t i = 0; i < a.size(); i++) res[i + 1] = (ll)a[i] * inv(i + 1) % MOD;
        return res;
    }

    /// q = P' / P from q[k] = P'[k] - (P[1] q[k - 1] + ... + P[k] q[0]), then log = integral of q
    static vector<int> log(vector<int> a, int n){
        if (n == 0) return {};
        a = cut(a, n);
        vector<int> d = cut(derivative(a), n - 1), q(n - 1, 0);
        for (int k = 0; k < n - 1; k++){
            ll s = d[k];
            for (int j = 1; j <= k; j++) s = (s - (ll)a[j] * q[k - j]) % MOD;
            q[k] = (s + MOD) % MOD;
        }
        return integral(q);
    }

    /// E' = P' E gives k E[k] = sum j P[j] E[k - j]
    static vector<int> exp(vector<int> a, int n){
        if (n == 0) return {};
        a = cut(a, n);
        vector<int> e(n, 0);
        e[0] = 1;
        for (int k = 1; k < n; k++){
            ll s = 0;
            for (int j = 1; j <= k; j++) s = (s + (ll)j * a[j] % MOD * e[k - j]) % MOD;
            e[k] = s * inv(k) % MOD;
        }
        return e;
    }

    /// R[0] = root, 2 R[0] R[k] = Y[k] - (R[1] R[k - 1] + ... + R[k - 1] R[1])
    static vector<int> sqrt_series(const vector<int>& a, int n, int i, ll root){
        vector<int> res(n, 0), y = cut(vector<int>(a.begin() + i, a.end()), n - i), r(n - i / 2, 0);
        y.resize(r.size(), 0);
        ll den = inv(2 * root % MOD);
        r[0] = root;
        for (size_t k = 1; k < r.size(); k++){
            ll s = y[k];
            for (size_t j = 1; j < k; j++) s = (s - (ll)r[j] * r[k - j]) % MOD;
            r[k] = (s + MOD) % MOD * den % MOD;
        }
        for (size_t k = 0; k < r.size(); k++) res[i / 2 + k] = r[k];
        return res;
    }

    static vector<int> multiply_cut(const vector<int>& x, const vector<int>& y, int n){
        vector<int> res(n, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; i + j < n; j++) res[i + j] = (res[i + j] + (ll)x[i] * y[j]) % MOD;
        }
        return res;
    }

    static vector<int> pow(vector<int> a, ll k, int n){
        vector<int> res(n, 0), base = cut(a, n);
        if (n == 0) return res;
        res[0] = 1;
        for (; k > 0; k >>= 1, base = multiply_cut(base, base, n)){
            if (k & 1) res = multiply_cut(res, base, n);
        }
        return res;
    }

    static int evaluate(const vector<int>& a, ll x){
        ll res = 0;
        for (int i = (int)a.size() - 1; i >= 0; i--) res = (res * x + a[i]) % MOD;
        return res;
    }

    /// Lagrange: sum y[i] * prod (x - x[j]) / (x[i] - x[j]), with M / (x - x[i]) by synthetic division
    static vector<int> interpolate(const vector<int>& xs, const vector<int>& ys){
        int n = xs.size();
        vector<int> m = {1}, res(n, 0);
        for (int i = 0; i < n; i++) m = multiply(m, {(MOD - xs[i]) % MOD, 1});

        for (int i = 0; i < n; i++){
            vector<int> q(n, 0);
            ll carry = 0, den = 1;
            for (int d = n; d >= 1; d--){
                carry = (m[d] + carry * xs[i]) % MOD;
                q[d - 1] = carry;
            }
            for (int j = 0; j < n; j++) if (j != i) den = den * ((xs[i] - xs[j] + MOD) % MOD) % MOD;

            ll c = ys[i] * inv(den) % MOD;
            for (int d = 0; d < n; d++) res[d] = (res[d] + c * q[d]) % MOD;
        }
        return trim(res);
    }

    /// Horner with (x + c) in place of x
    static vector<int> taylor_shift(const vector<int>& a, ll c){
        vector<int> res;
        for (int i = (int)a.size() - 1; i >= 0; i--){
            res = multiply(res, {(int)c, 1});
            if (res.empty()) res = {0};
            res[0] = (res[0] + a[i]) % MOD;
        }
        return cut(res, a.size());
    }
};

template<int MOD>
vector<int> random_coefficients(int n){
    vector<int> v(n);
    int style = stress::rand_int(0, 3);
    for (auto& x : v){
        if (style == 0) x = stress::rand_int(0, 2);
        else if (style == 1) x = stress::rand_int(0, 3) ? 0 : stress::rand_int(0, MOD - 1);
        else x = stress::rand_int(0, MOD - 1);
    }
    return v;
}

template<int MOD>
vector<int> distinct_points(int n){
    set<int> seen;
    vector<int> pts;
    while ((int)pts.size() < n){
        int x = stress::rand_int(0, MOD - 1);
        if (seen.insert(x).second) pts.push_back(x);
    }
    return pts;
}

/// Every operation on sizes n and m against the schoolbook versions, series lengths stay within limit = 2^(v - 1)
template<int MOD>
void check_all(int n, int m, int limit = 1 << 20){
    using P = Polynomial<MOD>;
    using S = Schoolbook<MOD>;

    vector<int> x = random_coefficients<MOD>(n), y = random_coefficients<MOD>(m);
    P px(x), py(y);

    /// Constructor reduction of negative and oversized values
    vector<ll> raw(n);
    for (int i = 0; i < n; i++) raw[i] = x[i] - (ll)MOD * stress::rand_int(-1000000, 1000000);
    assert(P(raw).a == x);

    assert((px * py).a == S::trim(S::multiply(x, y)));
    assert((px + py).a == S::add(x, y, 1));
    assert((px - py).a == S::add(x, y, -1));
    assert(px.derivative().a == S::derivative(x));
    assert(px.integral().a == S::integral(x));
    int c = stress::rand_int(0, MOD - 1);
    assert(px.taylor_shift(c).a == S::taylor_shift(x, c));
    assert(px.taylor_shift(c - (ll)MOD * stress::rand_int(1, 1000000)).a == S::taylor_shift(x, c));
    assert(px.evaluate(c - (ll)MOD * stress::rand_int(-1000000, 1000000)) == S::evaluate(x, c));

    if (!S::trim(y).empty()){
        auto [q, r] = px.divmod(py);
        auto [bq, br] = S::divmod(x, y);
        assert(q.a == bq && r.a == br && (px / py).a == bq && (px % py).a == br);
    }

    /// Series operations read only the first k coefficients, so x keeps its random tail past k
    int k = n == limit ? limit : min<int>(limit, stress::rand_int(0, n + 3));
    vector<int> unit = x;
    unit.resize(max<int>(n, 1));
    if (unit[0] == 0) unit[0] = stress::rand_int(1, MOD - 1);
    assert(P(unit).inverse(k).a == S::inverse(unit, k));

    unit[0] = 1;
    assert(P(unit).log(k).a == S::log(unit, k));

    unit[0] = 0;
    assert(P(unit).exp(k).a == S::exp(unit, k));

    /// pow: leading zeros, k both small (exhaustive) and huge
    vector<int> shifted = x;
    shifted.insert(shifted.begin(), stress::rand_int(0, 4), 0);
    for (ll e = 0; e <= 3; e++) assert(P(shifted).pow(e, k).a == S::pow(shifted, e, k));
    ll big = stress::rand_int(0, 1) && k <= 300 ? stress::rand_int(0, 1000000000000000000LL) : stress::rand_int(4, 40);
    assert(P(shifted).pow(big, k).a == S::pow(shifted, big, k));

    /// sqrt: half the time a square's lowest coefficient (root known), otherwise random
    vector<int> sq = shifted;
    int lead = 0, lim = min<int>(k, sq.size());
    while (lead < lim && sq[lead] == 0) lead++;
    ll known = -1;
    if (lead < lim && stress::rand_int(0, 1)){
        known = stress::rand_int(1, MOD - 1);
        sq[lead] = known * known % MOD;
        known = min<ll>(known, MOD - known);
    }

    auto root = P(sq).sqrt(k);
    if (lead == lim) assert(root && root->a == vector<int>(k, 0));
    else if (lead % 2 || S::pw(sq[lead], (MOD - 1) / 2) != 1) assert(!root);
    else{
        assert(root);
        ll s = root->a[lead / 2];
        assert(s != 0 && s * s % MOD == sq[lead] && s <= MOD - s && (known < 0 || s == known));
        assert(root->a == S::sqrt_series(S::cut(sq, k), k, lead, s));
    }

    /// Multipoint evaluation at arbitrary points (repeats allowed) and interpolation at distinct ones
    vector<int> pts(m);
    for (auto& p : pts) p = stress::rand_int(0, 2) ? stress::rand_int(0, MOD - 1) : stress::rand_int(0, 3);
    auto values = px.evaluate(pts);
    for (int i = 0; i < m; i++) assert(values[i] == S::evaluate(x, pts[i]));

    int cnt = min(m, MOD);
    auto xs = distinct_points<MOD>(cnt), ys = random_coefficients<MOD>(cnt);
    auto poly = P::interpolate(xs, ys);
    assert(poly.a == S::interpolate(xs, ys));
}

int main(){
    for (ll it = 0; it < stress::scaled(120); it++){
        int n = it < 40 ? it : stress::rand_int(0, it % 10 ? 80 : 300), m = stress::rand_int(0, it % 10 ? 80 : 300);
        check_all<998244353>(n, m);
    }

    for (ll it = 0; it < stress::scaled(25); it++){
        int n = stress::rand_int(0, 120), m = stress::rand_int(0, 120);
        check_all<167772161>(n, m);
        check_all<469762049>(n, m);
        check_all<754974721>(n, m);
    }

    /// The size limit 2^(v - 1) with v = 4 for 17, 5 for 97 and 12 for 12289: 17 and 97 push the series operations
    /// to that limit, 12289 with n = m = 2048 reaches the transform length 2^v through multiply and interpolate
    for (ll it = 0; it < stress::scaled(200); it++){
        check_all<17>(stress::rand_int(0, 8), stress::rand_int(0, 8), 8);
        check_all<97>(stress::rand_int(0, 16), stress::rand_int(0, 16), 16);
    }
    for (ll it = 0; it < stress::scaled(2); it++){
        check_all<12289>(2048, 2048, 2048);
    }

    return 0;
}
