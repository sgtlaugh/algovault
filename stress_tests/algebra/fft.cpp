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

/// Resident memory in kB, -1 if /proc is unavailable
long long resident_kb(){
    ifstream status("/proc/self/status");
    string line;
    while (getline(status, line)){
        if (line.rfind("VmRSS:", 0) == 0) return atoll(line.c_str() + 6);
    }
    return -1;
}

string random_bits(int n){
    string s(n, '0');
    for (auto& c : s) c = '0' + stress::rand_int(0, 1);
    return s;
}

int main(){
    using namespace fft;

    /// The static arrays span ~384 MB and must stay untouched until used, a non-constexpr
    /// ComplexNum constructor would dynamically initialise all of them at startup
    long long startup_kb = resident_kb();
    assert(startup_kb < 128 * 1024);

    /// Coefficients past LL_MULTIPLY_LIMIT^2 ~ 2.25e18 but below 2^63 must not wrap
    const ll big = LL_MULTIPLY_LIMIT - 1;
    assert(ll_multiply({big, big, big, big}, {big, big, big, big}) == naive({big, big, big, big}, {big, big, big, big}));
    assert(ll_multiply({big, big, big, big}, {big, big, big, big - 1}) == naive({big, big, big, big}, {big, big, big, big - 1}));

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, it % 20 ? 300 : 1200), m = stress::rand_int(1, it % 20 ? 300 : 1200);

        auto a = random_vector(n, -10000, 10000), b = random_vector(m, -10000, 10000);
        assert(multiply(a, b) == naive(a, b));
        assert(square(a) == naive(a, a));

        ll mod = stress::rand_int(2, 2147483647);
        auto p = random_vector(n, 0, mod - 1), q = random_vector(m, 0, mod - 1);
        assert(mod_multiply(p, q, mod) == naive(p, q, mod));
        assert(mod_multiply(p, p, mod) == naive(p, p, mod));  /// equal inputs skip the second transform

        /// ll_multiply is exact while every output coefficient stays below 2^63, hi is pushed to that edge half the time
        int xn = min(n, 1000), yn = min(m, 1000);
        ll hi = stress::rand_int(0, 1) ? 30000000 : min<ll>(LL_MULTIPLY_LIMIT - 1, sqrtl(9.2e18L / max(xn, yn)));
        auto x = random_vector(xn, 0, hi), y = random_vector(yn, 0, hi);
        assert(ll_multiply(x, y) == naive(x, y));
        assert(ll_multiply(x, x) == naive(x, x));

        int k = min(n, 500);
        auto c1 = random_vector(k, -1000, 1000), c2 = random_vector(k, -1000, 1000);
        assert(convolution(c1, c2) == naive_circular(c1, c2));
        auto d1 = random_vector(k, 0, mod - 1), d2 = random_vector(k, 0, mod - 1);
        assert(mod_convolution(d1, d2, mod) == naive_circular(d1, d2, mod));
        auto e1 = random_vector(k, 0, 1000000), e2 = random_vector(k, 0, 1000000);
        assert(ll_convolution(e1, e2) == naive_circular(e1, e2));

        string text = random_bits(stress::rand_int(1, 400)), pattern = random_bits(stress::rand_int(1, text.size()));
        auto hd = hamming_distance(text.c_str(), pattern.c_str()), ac = and_convolution(text.c_str(), pattern.c_str());
        assert(hd.size() == text.size() - pattern.size() + 1 && ac.size() == hd.size());
        for (size_t i = 0; i < hd.size(); i++){
            ll mismatches = 0, common = 0;
            for (size_t j = 0; j < pattern.size(); j++) mismatches += text[i + j] != pattern[j], common += text[i + j] == '1' && pattern[j] == '1';
            assert(hd[i] == mismatches && ac[i] == common);
        }
    }

    /// Near the precision limit: n = 10^4 with values up to 2 * 10^6 is exact with a long double pi (limit ~1.6 * 10^7)
    /// but not with a double pi (limit ~7.5 * 10^5), 200 sampled coefficients are checked exactly
    for (long long it = 0; it < stress::scaled(2); it++){
        int n = 10000;
        auto a = random_vector(n, 0, 2000000), b = random_vector(n, 0, 2000000);
        auto c = multiply(a, b);
        for (int t = 0; t < 200; t++){
            int k = stress::rand_int(0, 2 * n - 2);
            __int128 expected = 0;
            for (int i = max(0, k - n + 1); i <= min(k, n - 1); i++) expected += (__int128)a[i] * b[k - i];
            assert(c[k] == (ll)expected);
        }
    }

    return 0;
}
