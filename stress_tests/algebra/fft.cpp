#include "../common.h"

#define main library_main
#include "../../code_library/algebra/fft.cpp"
#undef main

typedef long long ll;

vector<ll> random_vector(int n, ll lo, ll hi){
    vector<ll> v(n);
    for (auto& x : v) x = stress::rand_int(lo, hi);
    return v;
}

/// Exact reference, __int128 keeps ll_multiply sized products exact before the optional reduction
vector<ll> naive(const vector<ll>& a, const vector<ll>& b, ll mod = 0){
    vector<__int128> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) c[i + j] += (__int128)a[i] * b[j];

    vector<ll> res(c.size());
    for (size_t i = 0; i < c.size(); i++) res[i] = mod ? (ll)(((c[i] % mod) + mod) % mod) : (ll)c[i];
    return res;
}

vector<ll> naive_circular(const vector<ll>& a, const vector<ll>& b, ll mod = 0){
    int n = a.size();
    vector<__int128> c(n, 0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) c[(i + j) % n] += (__int128)a[i] * b[j];

    vector<ll> res(n);
    for (int i = 0; i < n; i++) res[i] = mod ? (ll)(c[i] % mod) : (ll)c[i];
    return res;
}

/// The header's circular convolution recipe: fold the linear product back onto its first n coefficients
vector<ll> fold(vector<ll> c, int n, ll mod = 0){
    for (int i = n; i < (int)c.size(); i++) c[i - n] = mod ? (c[i - n] + c[i]) % mod : c[i - n] + c[i];
    c.resize(n);
    return c;
}

/// All-equal inputs have a closed form, c[k] = va * vb * (number of pairs i + j = k), so large n needs no brute force
bool equal_product_is_exact(FFT& fft, int n, int m, ll va, ll vb){
    auto c = n == m && va == vb ? fft.square(vector<ll>(n, va)) : fft.multiply(vector<ll>(n, va), vector<ll>(m, vb));
    for (int k = 0; k < n + m - 1; k++){
        if (c[k] != va * vb * min({k + 1, n, m, n + m - 1 - k})) return false;
    }
    return true;
}

string random_bits(int n){
    string s(n, '0');
    for (auto& c : s) c = '0' + stress::rand_int(0, 1);
    return s;
}

int main(){
    FFT fft;

    /// Coefficients past LL_MULTIPLY_LIMIT^2 ~ 2.25e18 but below 2^63 must not wrap
    const ll big = FFT::LL_MULTIPLY_LIMIT - 1;
    assert(fft.ll_multiply({big, big, big, big}, {big, big, big, big}) == naive({big, big, big, big}, {big, big, big, big}));
    assert(fft.ll_multiply({big, big, big, big}, {big, big, big, big - 1}) == naive({big, big, big, big}, {big, big, big, big - 1}));

    assert(fft.multiply({}, {}).empty() && fft.square({}).empty() && fft.ll_multiply({1}, {}).empty() && fft.mod_multiply({}, {3}, 5).empty());

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, it % 20 ? 300 : 1200), m = stress::rand_int(1, it % 20 ? 300 : 1200);

        auto a = random_vector(n, -10000, 10000), b = random_vector(m, -10000, 10000);
        assert(fft.multiply(a, b) == naive(a, b));
        assert(fft.square(a) == naive(a, a));

        /// Past EXACT_LIMIT multiply and square switch to the split product, signed values up to LL_MULTIPLY_LIMIT,
        /// hi keeps every coefficient below 2^63
        ll big_hi = FFT::LL_MULTIPLY_LIMIT - 1, small_hi = 9e18L / max(n, m) / big_hi;
        auto sa = random_vector(n, -big_hi, big_hi), sb = random_vector(m, -small_hi, small_hi);
        assert(fft.multiply(sa, sb) == naive(sa, sb));
        assert(fft.multiply(sb, sa) == naive(sb, sa));
        ll square_hi = min<ll>(big_hi, sqrtl(9e18L / m));
        auto sq = random_vector(m, it % 3 ? -square_hi : 0, it % 3 == 1 ? 0 : square_hi);
        assert(fft.square(sq) == naive(sq, sq));
        assert(fft.multiply(sq, sq) == naive(sq, sq));

        ll mod = stress::rand_int(2, 2147483647);
        auto p = random_vector(n, it % 4 ? 0 : -5 * mod, it % 4 ? mod - 1 : 5 * mod), q = random_vector(m, 0, mod - 1);
        assert(fft.mod_multiply(p, q, mod) == naive(p, q, mod));
        assert(fft.mod_multiply(p, p, mod) == naive(p, p, mod));  /// equal inputs skip the second transform

        /// ll_multiply is exact while every output coefficient stays below 2^63, hi is pushed to that edge half the time
        int xn = min(n, 1000), yn = min(m, 1000);
        ll hi = stress::rand_int(0, 1) ? 30000000 : min<ll>(FFT::LL_MULTIPLY_LIMIT - 1, sqrtl(9.2e18L / max(xn, yn)));
        auto x = random_vector(xn, 0, hi), y = random_vector(yn, 0, hi);
        assert(fft.ll_multiply(x, y) == naive(x, y));
        assert(fft.ll_multiply(x, x) == naive(x, x));

        int k = min(n, 500);
        auto c1 = random_vector(k, -1000, 1000), c2 = random_vector(k, -1000, 1000);
        assert(fold(fft.multiply(c1, c2), k) == naive_circular(c1, c2));
        auto d1 = random_vector(k, 0, mod - 1), d2 = random_vector(k, 0, mod - 1);
        assert(fold(fft.mod_multiply(d1, d2, mod), k, mod) == naive_circular(d1, d2, mod));
        auto e1 = random_vector(k, 0, 1000000), e2 = random_vector(k, 0, 1000000);
        assert(fold(fft.ll_multiply(e1, e2), k) == naive_circular(e1, e2));

        string text = random_bits(stress::rand_int(0, 400)), pattern = random_bits(stress::rand_int(0, 9) ? stress::rand_int(0, text.size()) : 0);
        auto hd = fft.hamming_distance(text, pattern), ac = fft.common_ones(text, pattern);
        assert(hd.size() == text.size() - pattern.size() + 1 && ac.size() == hd.size());
        for (size_t i = 0; i < hd.size(); i++){
            ll mismatches = 0, common = 0;
            for (size_t j = 0; j < pattern.size(); j++) mismatches += text[i + j] != pattern[j], common += text[i + j] == '1' && pattern[j] == '1';
            assert(hd[i] == mismatches && ac[i] == common);
        }
    }

    /// One transform was wrong on each of these: the error follows n_max * max|value|^2, not the output bound
    /// n_min * max|a| * max|b|, which is only 1.8e13 for the first one
    assert(equal_product_is_exact(fft, 30000, 30000, 2, 300000000));
    assert(equal_product_is_exact(fft, 30000, 30000, -2, -300000000));
    assert(equal_product_is_exact(fft, 1000, 30000, 300000000, -1));
    assert(equal_product_is_exact(fft, 20000, 20000, 7071068, 7071068));
    assert(equal_product_is_exact(fft, 20000, 20000, -7071068, 7071067));
    assert(equal_product_is_exact(fft, 1, 100000, 1, 100000000));
    assert(equal_product_is_exact(fft, 100000, 1, -100000000, 1));

    /// Both sides of each switch: double up to DOUBLE_LIMIT, long double up to EXACT_LIMIT, the split product past it,
    /// at the largest all-equal value on the cheaper side and one more; also the length 2 transform with a half of 1
    for (int n : {100000, 1000000}){
        for (long double limit : {(long double)FFT::DOUBLE_LIMIT, FFT::EXACT_LIMIT}){
            ll v = sqrtl(limit / n);
            while ((long double)n * (v + 1) * (v + 1) <= limit) v++;
            while ((long double)n * v * v > limit) v--;
            assert(equal_product_is_exact(fft, n, n, v, v));
            assert(equal_product_is_exact(fft, n, n, v + 1, -(v + 1)));
            assert(equal_product_is_exact(fft, n, n / 3, v, v + 1));
        }
    }
    /// Between double's first error (about 3e14) and EXACT_LIMIT only long double is exact
    assert(equal_product_is_exact(fft, 100000, 100000, 100000, 100000));
    assert(equal_product_is_exact(fft, 1000000, 1000000, 40000, -40000));
    assert(fft.multiply({7}, {-6}) == vector<ll>({-42}));
    assert(fft.square({-9}) == vector<ll>({81}));

    /// Near the precision limit: n = 10^4 with values up to 2 * 10^6 is exact with a long double pi (limit ~1.6 * 10^7)
    /// but not with a double pi (limit ~7.5 * 10^5), 200 sampled coefficients are checked exactly
    for (long long it = 0; it < stress::scaled(2); it++){
        int n = 10000;
        auto a = random_vector(n, 0, 2000000), b = random_vector(n, 0, 2000000);
        auto c = fft.multiply(a, b);
        for (int t = 0; t < 200; t++){
            int k = stress::rand_int(0, 2 * n - 2);
            __int128 expected = 0;
            for (int i = max(0, k - n + 1); i <= min(k, n - 1); i++) expected += (__int128)a[i] * b[k - i];
            assert(c[k] == (ll)expected);
        }
    }

    return 0;
}
