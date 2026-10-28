#include "../common.h"

#define main library_main
#include "../../code_library/algebra/online_convolution.cpp"
#undef main

#define main ntt_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

using Next = function<long long(int, long long)>;
using Multiply = function<vector<long long>(const vector<long long>&, const vector<long long>&)>;

long long power_mod(long long b, long long e, long long m){
    long long res = 1 % m;
    for (b %= m; e; e >>= 1, b = b * b % m){
        if (e & 1) res = res * b % m;
    }

    return res;
}

/// The recurrence term by term in O(n^2), reading g_1..g_i only
vector<long long> brute(int n, vector<long long> g, long long mod, const Next& next){
    g.resize(n, 0);
    for (auto& x : g) x = (x % mod + mod) % mod;

    vector<long long> f(n);
    for (int i = 0; i < n; i++){
        long long s = 0;
        for (int j = 0; j < i; j++) s = (s + f[j] * g[i - j]) % mod;
        f[i] = (next(i, s) % mod + mod) % mod;
    }

    return f;
}

/// Linear, nonlinear with negative values, Bell's i b_i = s_i (prime mod only), and one ignoring s
Next make_next(int kind, long long mod, int n){
    long long a = stress::rand_int(0, mod - 1);
    vector<long long> noise(n);
    for (auto& x : noise) x = stress::rand_int(-(1LL << 40), 1LL << 40);

    if (kind == 0) return [=](int i, long long s){ return (s * a + i) % mod; };
    if (kind == 1) return [=](int i, long long s){ return s * s % mod - noise[i]; };
    if (kind == 2) return [=](int i, long long s){ return i ? s * power_mod(i, mod - 2, mod) % mod : 1; };
    return [=](int i, long long){ return noise[i]; };
}

vector<long long> random_g(int n){
    vector<long long> g(stress::rand_int(0, n + 3));
    for (auto& x : g){
        int kind = stress::rand_int(0, 9);
        x = kind == 0 ? LLONG_MIN : kind == 1 ? LLONG_MAX : stress::rand_int(-(1LL << 50), 1LL << 50);
    }

    return g;
}

bool is_prime(long long x){
    if (x < 2) return false;
    for (long long d = 2; d * d <= x; d++){
        if (x % d == 0) return false;
    }

    return true;
}

void check(int n, long long mod, const Multiply& multiply){
    int kind = stress::rand_int(0, 3);
    if (kind == 2 && !(is_prime(mod) && n < mod)) kind = 0;
    Next next = make_next(kind, mod, n);
    vector<long long> g = random_g(n);

    /// next must see i = 0, 1, ..., n - 1 once each, so stateful callers work
    int calls = 0;
    auto counted = [&](int i, long long s){
        assert(i == calls++);
        assert(0 <= s && s < mod);
        return next(i, s);
    };

    assert(online_convolution(n, g, mod, multiply, counted) == brute(n, g, mod, next));
    assert(calls == n);
}

/// f_i = next(i, s_i) at scale, with s read off one full product f * g where g_0 = 0
void check_large(int n, long long mod, const Multiply& multiply){
    Next next = make_next(stress::rand_int(0, 1), mod, n);
    vector<long long> g(n);
    for (auto& x : g) x = stress::rand_int(0, mod - 1);

    vector<long long> f = online_convolution(n, g, mod, multiply, next);
    g[0] = 0;
    vector<long long> s = multiply(f, g);
    for (int i = 0; i < n; i++) assert(f[i] == (next(i, s[i]) % mod + mod) % mod);
}

int main(){
    using Poly = vector<long long>;
    auto ntt998 = [](const Poly& a, const Poly& b){ return ntt::multiply(a, b); };
    auto ntt7340033 = [](const Poly& a, const Poly& b){ return ntt::multiply<7340033>(a, b); };
    auto mod_mul = [](long long m){
        return Multiply([m](const Poly& a, const Poly& b){ return ntt::mod_multiply(a, b, m); });
    };
    /// Schoolbook product with trailing zeros trimmed, like Polynomial's operator*, so prod can be short
    auto trimmed = [](long long m){
        return Multiply([m](const Poly& a, const Poly& b){
            Poly res(a.size() + b.size() - 1, 0);
            for (size_t i = 0; i < a.size(); i++){
                for (size_t j = 0; j < b.size(); j++) res[i + j] = (res[i + j] + a[i] * b[j]) % m;
            }
            while (!res.empty() && res.back() == 0) res.pop_back();
            return res;
        });
    };

    /// Every n across the leaf size and two levels above it, for NTT primes, a CRT prime, composites and mod = 1
    vector<long long> moduli = {1, 2, 6, 1000000000, 1000000007, (1LL << 31) - 1};
    for (int n = 0; n <= 100; n++){
        check(n, 998244353, ntt998);
        check(n, 7340033, ntt7340033);
        check(n, moduli[n % moduli.size()], mod_mul(moduli[n % moduli.size()]));
        check(n, 5, trimmed(5));
    }

    /// Random sizes and moduli, the zero-heavy mod 2 and 3 cases make trimmed products short
    for (long long it = 0; it < stress::scaled(15); it++){
        int n = stress::rand_int(0, 1500);
        long long m = stress::rand_int(0, 2) ? stress::rand_int(1, (1LL << 31) - 1) : stress::rand_int(1, 50);
        check(n, 998244353, ntt998);
        check(n, m, mod_mul(m));

        long long small = stress::rand_int(2, 3);
        check(stress::rand_int(0, 300), small, trimmed(small));
    }

    check_large(1 << 16, 998244353, ntt998);
    check_large(30000, (1LL << 31) - 1, mod_mul((1LL << 31) - 1));

    return 0;
}
