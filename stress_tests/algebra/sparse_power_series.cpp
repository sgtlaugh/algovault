#include "../common.h"

#define main library_main
#include "../../code_library/algebra/sparse_power_series.cpp"
#undef main
#define main polynomial_main
#include "../../code_library/algebra/polynomial.cpp"
#undef main

using Series = vector<long long>;

long long power_mod(long long b, long long e, long long mod){
    long long res = 1;
    b %= mod;
    while (e){
        if (e & 1) res = res * b % mod;
        b = b * b % mod;
        e >>= 1;
    }

    return res;
}

Series multiply(const Series& a, const Series& b, int n, long long mod){
    Series c(n, 0);
    for (int i = 0; i < n && i < (int)a.size(); i++){
        if (!a[i]) continue;
        for (int j = 0; i + j < n && j < (int)b.size(); j++) c[i + j] = (c[i + j] + a[i] * b[j]) % mod;
    }

    return c;
}

Series dense(const vector<pair<int, long long>>& raw, int n, long long mod){
    Series f(n, 0);
    for (auto [d, c] : raw){
        if (d < n) f[d] = ((f[d] + c) % mod + mod) % mod;
    }

    return f;
}

/// Square and multiply on dense series, so huge exponents and leading zeros come out without any shortcut
Series brute_pow(const Series& f, long long k, int n, long long mod){
    Series res(n, 0), b = f;
    if (n) res[0] = 1;
    while (k){
        if (k & 1) res = multiply(res, b, n, mod);
        b = multiply(b, b, n, mod);
        k >>= 1;
    }

    return res;
}

/// Long division with the dense series, independent of the sparse recurrence
Series brute_inverse(const Series& f, int n, long long mod){
    Series g(n, 0);
    long long c0_inv = power_mod(f[0], mod - 2, mod);
    for (int i = 0; i < n; i++){
        long long rest = i == 0 ? 1 : 0;
        for (int j = 1; j <= i; j++) rest = (rest - f[j] * g[i - j]) % mod;
        g[i] = (rest + mod) % mod * c0_inv % mod;
    }

    return g;
}

/// e^f = sum f^j / j!, and ln f = sum (-1)^(j + 1) (f - 1)^j / j: both stop at j = n since f - 1 and f start at x
Series brute_exp(const Series& f, int n, long long mod){
    Series res(n, 0), term(n, 0);
    if (!n) return res;
    term[0] = 1;
    for (int j = 0; j < n; j++){
        for (int i = 0; i < n; i++) res[i] = (res[i] + term[i]) % mod;
        term = multiply(term, f, n, mod);
        long long inv = power_mod(j + 1, mod - 2, mod);
        for (auto& t : term) t = t * inv % mod;
    }

    return res;
}

Series brute_log(const Series& f, int n, long long mod){
    Series res(n, 0), u = f, term(n, 0);
    if (!n) return res;
    u[0] = 0;
    term[0] = 1;
    for (int j = 1; j < n; j++){
        term = multiply(term, u, n, mod);
        long long w = power_mod(j, mod - 2, mod);
        if (j % 2 == 0) w = (mod - w) % mod;
        for (int i = 0; i < n; i++) res[i] = (res[i] + term[i] * w) % mod;
    }

    return res;
}

bool is_residue(long long x, long long mod){
    return power_mod(x, (mod - 1) / 2, mod) == 1;
}

vector<int> as_int(const Series& v){
    return vector<int>(v.begin(), v.end());
}

/// A sparse series: few terms at random degrees, sometimes repeated, some coefficients negative or multiples of mod
vector<pair<int, long long>> random_terms(int max_degree, int k, long long mod){
    vector<pair<int, long long>> raw;
    for (int j = 0; j < k; j++){
        int d = stress::rand_int(0, max_degree);
        long long c = stress::rand_int(0, 3) ? stress::rand_int(-3 * mod, 3 * mod) : mod * stress::rand_int(-2, 2);
        raw.push_back({d, c});
    }

    return raw;
}

/// exp needs a zero constant term and log a constant term 1, so force them on a copy
void check_exp_log(const vector<pair<int, long long>>& raw, int n, int mod){
    auto raw0 = raw;
    raw0.erase(remove_if(raw0.begin(), raw0.end(), [](auto& t){ return t.first == 0; }), raw0.end());
    SparseSeries e(raw0, mod);
    Series ed = dense(raw0, n, mod);
    assert(e.exp(n) == as_int(brute_exp(ed, n, mod)));

    auto raw1 = raw0;
    raw1.push_back({0, 1});
    SparseSeries l(raw1, mod);
    Series ld = dense(raw1, n, mod);
    assert(l.log(n) == as_int(brute_log(ld, n, mod)));
}

void check_small(int mod, int max_n = 24){
    int n = stress::rand_int(0, min(mod, max_n));
    int k = stress::rand_int(0, 6);
    auto raw = random_terms(n + 3, k, mod);
    SparseSeries f(raw, mod);
    Series fd = dense(raw, n, mod);

    for (int j = 1; j < (int)f.terms.size(); j++) assert(f.terms[j - 1].first < f.terms[j].first);
    for (auto [d, c] : f.terms) assert(1 <= c && c < mod && d >= 0);

    /// Terms of degree >= n must not change any result
    vector<pair<int, long long>> low;
    for (auto [d, c] : raw) if (d < n) low.push_back({d, c});
    SparseSeries g(low, mod);

    long long k_exp = stress::rand_int(0, 3) ? stress::rand_int(0, 12) : stress::rand_int(0, (long long)4e18);
    Series want = brute_pow(fd, k_exp, n, mod);
    assert(f.pow(k_exp, n) == as_int(want));
    assert(g.pow(k_exp, n) == as_int(want));

    if (n > 0 && fd[0]){
        assert(f.inverse(n) == as_int(brute_inverse(fd, n, mod)));
    }

    check_exp_log(raw, n, mod);
    if (mod == 2) return;

    auto root = f.sqrt(n);
    int s = 0;
    while (s < n && fd[s] == 0) s++;
    bool exists = s == n || (s % 2 == 0 && is_residue(fd[s], mod));
    assert(root.has_value() == exists);
    if (root){
        Series r(root->begin(), root->end());
        assert(multiply(r, r, n, mod) == fd);
        if (s < n) assert((*root)[s / 2] <= mod - (*root)[s / 2]);
        assert(g.sqrt(n) == root);
    }
}

/// Large n against the NTT Polynomial, a separate implementation of the same operations
void check_large(int n, int k){
    using P = Polynomial<998244353>;
    const int mod = 998244353;
    auto raw = random_terms(n - 1, k, mod);
    raw.push_back({0, stress::rand_int(1, mod - 1)});
    SparseSeries f(raw, mod);
    P p(dense(raw, n, mod));

    assert(f.inverse(n) == p.inverse(n).a);
    long long k_exp = stress::rand_int(0, (long long)1e18);
    assert(f.pow(k_exp, n) == p.pow(k_exp, n).a);
    assert(f.sqrt(n) == (p.sqrt(n) ? optional<vector<int>>(p.sqrt(n)->a) : nullopt));

    auto raw0 = raw;
    raw0.erase(remove_if(raw0.begin(), raw0.end(), [](auto& t){ return t.first == 0; }), raw0.end());
    assert(SparseSeries(raw0, mod).exp(n) == P(dense(raw0, n, mod)).exp(n).a);

    raw0.push_back({0, 1});
    assert(SparseSeries(raw0, mod).log(n) == P(dense(raw0, n, mod)).log(n).a);

    /// A shifted series: lowest degree s > 0 exercises the pow shift and the sqrt parity rule
    auto shifted = raw;
    int s = stress::rand_int(1, 40);
    for (auto& t : shifted) t.first += s;
    SparseSeries fs(shifted, mod);
    P ps(dense(shifted, n, mod));
    long long small_k = stress::rand_int(0, n / s + 2);
    assert(fs.pow(small_k, n) == ps.pow(small_k, n).a);
    assert(fs.sqrt(n) == (ps.sqrt(n) ? optional<vector<int>>(ps.sqrt(n)->a) : nullopt));
}

int main(){
    const int primes[] = {3, 5, 7, 13, 101, 65537, 998244353, 1000000007, 2147483647};

    /// n = 0 and the zero series: every call returns n zeros, 0^0 = 1, and sqrt(0) = 0
    for (int mod : primes){
        SparseSeries zero({}, mod);
        assert(zero.pow(0, 3) == (vector<int>{1, 0, 0}));
        assert(zero.pow(5, 3) == (vector<int>{0, 0, 0}));
        assert(zero.sqrt(3) == (vector<int>{0, 0, 0}));
        assert(zero.exp(2) == (vector<int>{1, 0}));
        assert(zero.exp(0).empty());
        SparseSeries one({{0, 1}}, mod);
        assert(one.inverse(0).empty() && one.log(0).empty() && one.pow(7, 0).empty());
        assert(one.log(3) == (vector<int>{0, 0, 0}));
    }

    /// n = mod is the largest allowed: inverses of 1..mod - 1 all exist
    for (int mod : {3, 5, 7}){
        SparseSeries x({{1, 1}}, mod);
        assert((int)x.exp(mod).size() == mod);
    }

    for (long long it = 0; it < stress::scaled(4000); it++) check_small(primes[stress::rand_int(0, 8)]);

    /// mod = 2 allows n <= 2 and every call but sqrt; primes near 2^31 at n in the hundreds stress the 64-bit products
    for (long long it = 0; it < stress::scaled(200); it++) check_small(2);
    for (long long it = 0; it < stress::scaled(4); it++) check_small(stress::rand_int(0, 1) ? 2147483647 : 2147483629, 300);

    for (long long it = 0; it < stress::scaled(30); it++) check_large(stress::rand_int(1, 3000), stress::rand_int(1, 10));
    check_large(200000, 10);

    /// A broken precondition aborts in the library instead of returning garbage
    stress::expect_assert("algebra/sparse_power_series.cpp", []{ SparseSeries({{0, 2}}).log(4); });
    stress::expect_assert("algebra/sparse_power_series.cpp", []{ SparseSeries({{0, 1}}).exp(4); });
    stress::expect_assert("algebra/sparse_power_series.cpp", []{ SparseSeries({{1, 1}}).inverse(4); });
    stress::expect_assert("algebra/sparse_power_series.cpp", []{ SparseSeries({{1, 1}}, 7).exp(8); });
    stress::expect_assert("algebra/sparse_power_series.cpp", []{ SparseSeries({{-1, 1}}); });

    return 0;
}
