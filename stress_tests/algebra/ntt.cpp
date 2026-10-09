#include "../common.h"

#define main library_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

typedef long long ll;

vector<ll> random_vector(int n, ll lo, ll hi){
    vector<ll> v(n);
    for (auto& x : v) x = stress::rand_int(lo, hi);
    return v;
}

/// Exact O(nm) product in __int128, reduced into [0, mod) when mod is given
vector<ll> naive(const vector<ll>& a, const vector<ll>& b, ll mod = 0){
    if (a.empty() || b.empty()) return {};

    vector<__int128> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); i++){
        for (size_t j = 0; j < b.size(); j++){
            __int128 x = mod ? a[i] % mod : a[i], y = mod ? b[j] % mod : b[j];
            c[i + j] += x * y;
            if (mod) c[i + j] %= mod;
        }
    }

    vector<ll> res(c.size());
    for (size_t i = 0; i < c.size(); i++) res[i] = mod ? (ll)((c[i] % mod + mod) % mod) : (ll)c[i];
    return res;
}

/// Compares one prime's multiply against the naive product, sizes from 0 up to the given caps
template<unsigned MOD, unsigned ROOT>
void check_prime(int max_n, int max_m){
    int n = stress::rand_int(0, max_n), m = stress::rand_int(0, max_m);
    ll lo = stress::rand_int(0, 3) ? 0 : LLONG_MIN, hi = lo ? LLONG_MAX : MOD - 1;
    auto a = random_vector(n, lo, hi), b = random_vector(m, lo, hi);
    assert((ntt::multiply<MOD, ROOT>(a, b) == naive(a, b, MOD)));
}

/// Product of all-ones vectors of lengths n and m, c[k] = number of (i, j) with i + j = k
vector<ll> ones_product(int n, int m){
    vector<ll> c(n + m - 1);
    for (int k = 0; k < n + m - 1; k++) c[k] = min({k + 1, n, m, n + m - 1 - k});
    return c;
}

int main(){
    using namespace ntt;

    /// Exhaustive shapes on tiny inputs, every length pair up to 9 x 9 including empty
    for (int n = 0; n <= 9; n++){
        for (int m = 0; m <= 9; m++){
            for (int rep = 0; rep < 20; rep++){
                auto a = random_vector(n, -40, 40), b = random_vector(m, -40, 40);
                assert(multiply(a, b) == naive(a, b, 998244353));
                if (n + m - 1 <= 16) assert((multiply<17, 3>(a, b) == naive(a, b, 17)));
                assert(mod_multiply(a, b, rep + 1) == naive(a, b, rep + 1));
                assert(exact_multiply(a, b) == naive(a, b));
            }
        }
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        int big = it % 20 == 0;
        check_prime<998244353, 3>(big ? 3000 : 300, big ? 3000 : 300);
        check_prime<7340033, 3>(300, 300);
        check_prime<P1, 3>(300, 300);
        check_prime<P2, 3>(300, 300);
        check_prime<P3, 11>(300, 300);
        check_prime<257, 3>(128, 128);
        check_prime<17, 3>(8, 8);

        /// mod_multiply over any modulus, inputs anywhere in long long, the top modulus 2^31 - 1 a quarter of the time
        int n = stress::rand_int(0, big ? 2000 : 200), m = stress::rand_int(0, big ? 2000 : 200);
        ll mod = stress::rand_int(0, 3) ? stress::rand_int(1, (1LL << 31) - 1) : (1LL << 31) - 1;
        auto a = random_vector(n, LLONG_MIN, LLONG_MAX), b = random_vector(m, LLONG_MIN, LLONG_MAX);
        assert(mod_multiply(a, b, mod) == naive(a, b, mod));
        auto p = random_vector(n, mod - 2 > 0 ? mod - 2 : 0, mod - 1), q = random_vector(m, mod - 2 > 0 ? mod - 2 : 0, mod - 1);
        assert(mod_multiply(p, q, mod) == naive(p, q, mod));

        /// exact_multiply with signed values sized so the largest coefficient sits near 2^63 half the time
        int k = max(1, min(n, m));
        ll bound = stress::rand_int(0, 1) ? stress::rand_int(1, 1000000) : (ll)sqrtl(9.2e18L / k);
        auto x = random_vector(n, -bound, bound), y = random_vector(m, -bound, bound);
        assert(exact_multiply(x, y) == naive(x, y));
    }

    /// Length boundaries: exactly 2^8 for 257 and 2^20 for 7340033, all-ones products have a closed form
    auto u = random_vector(129, 0, 256), v = random_vector(128, 0, 256);
    assert((multiply<257, 3>(u, v) == naive(u, v, 257)));
    assert((multiply<7340033, 3>(vector<ll>((1 << 19) + 1, 1), vector<ll>(1 << 19, 1)) == ones_product((1 << 19) + 1, 1 << 19)));

    /// The header's maximum lengths, nightly only (STRESS_SCALE > 1) since a 2^24 three-prime product is slow under sanitizers
    /// (m - 1)^2 = 1 modulo m, so all-(m - 1) inputs reduce to the all-ones closed form while every true coefficient nears M
    if (stress::scaled(1) > 1){
        auto c = ones_product((1 << 22) + 1, 1 << 22);
        assert((multiply(vector<ll>((1 << 22) + 1, 1), vector<ll>(1 << 22, 1)) == c));

        const ll top = (1LL << 31) - 1;
        auto d = ones_product((1 << 23) + 1, 1 << 23);
        auto e = mod_multiply(vector<ll>((1 << 23) + 1, top - 1), vector<ll>(1 << 23, top - 1), top);
        for (auto& x : d) x %= top;
        assert(e == d);
    }

    /// Garner reconstruction up to M ~ 5.95e25, single products x * y < M checked exactly, x >= 6.5e6 keeps y in long long
    const unsigned __int128 M = (unsigned __int128)P1 * P2 * P3;
    for (long long it = 0; it < stress::scaled(2000); it++){
        ll x = stress::rand_int(6500000, 7713742000000LL), top = (ll)((M - 1) / x), y = stress::rand_int(top / 2, top);
        auto c = crt_multiply({x}, {y});
        assert(c.size() == 1 && c[0] == (unsigned __int128)x * y);
    }

    /// Exact boundary values of [-2^63, 2^63)
    assert((exact_multiply({3037000499LL, -3037000499LL}, {3037000499LL}) == vector<ll>{9223372030926249001LL, -9223372030926249001LL}));
    assert((exact_multiply({1LL << 62, 1LL << 62}, {1, -1}) == vector<ll>{1LL << 62, 0, -(1LL << 62)}));
    assert((exact_multiply({LLONG_MIN, LLONG_MAX, -1}, {1}) == vector<ll>{LLONG_MIN, LLONG_MAX, -1}));

    return 0;
}
