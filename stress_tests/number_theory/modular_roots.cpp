#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/modular_roots.cpp"
#undef main

long long brute_order(long long a, long long m){
    a %= m;
    if (gcd(a, m) != 1) return -1;
    long long x = 1, e = a % m;
    for (; e != 1 % m; x++) e = e * a % m;
    return x;
}

long long brute_primitive_root(long long m){
    long long phi = 0;
    for (long long x = 0; x < m; x++) phi += gcd(x, m) == 1;
    for (long long g = 0; g < m; g++){
        if (brute_order(g, m) == phi) return g;
    }
    return -1;
}

bool brute_is_prime(long long n){
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; d++){
        if (n % d == 0) return false;
    }
    return true;
}

/// x^k mod p for every x through discrete logs to a brute force generator, independent of pow_mod
vector<long long> brute_powers(long long k, long long p){
    long long g = brute_primitive_root(p);
    vector<long long> g_pow(p - 1, 1), res(p, k == 0);
    for (long long i = 1; i < p - 1; i++) g_pow[i] = g_pow[i - 1] * g % p;
    for (long long i = 0; i < p - 1; i++) res[g_pow[i]] = g_pow[k % (p - 1) * i % (p - 1)];
    return res;
}

void check_kth_roots(long long k, long long p){
    vector<long long> powers = brute_powers(k, p);
    vector<vector<long long>> roots_of(p);
    for (long long x = 0; x < p; x++) roots_of[powers[x]].push_back(x);

    for (long long a = 0; a < p; a++){
        const vector<long long>& expected = roots_of[a];
        long long x = kth_root(a, k, p);
        assert(expected.empty() ? x == -1 : 0 <= x && x < p && powers[x] == a);
        assert(kth_roots(a, k, p) == expected);
    }
}

long long random_prime(long long lo, long long hi){
    while (1){
        long long n = stress::rand_int(lo, hi);
        if (is_prime(n)) return n;
    }
}

/// Order, primitive root, square root and k-th roots against direct scans over every residue for small moduli,
/// then defining properties on random 62-bit primes
int main(){
    vector<long long> primes;
    for (long long n = 2; n < 1000; n++){
        if (brute_is_prime(n)) primes.push_back(n);
    }

    for (long long n = 1; n <= 30000; n++){
        long long prod = 1;
        for (auto [q, e] : factorize(n)){
            assert(brute_is_prime(q));
            for (int i = 0; i < e; i++) prod *= q;
        }
        assert(prod == n);
        assert(is_prime(n) == brute_is_prime(n));
    }
    for (long long it = 0; it < stress::scaled(300); it++){
        long long n = stress::rand_int(1, it % 2 ? 1000000000000LL : (1LL << 62) - 1), rest = n;
        for (auto [q, e] : factorize(n)){
            assert(is_prime(q));
            for (int i = 0; i < e; i++) assert(rest % q == 0), rest /= q;
        }
        assert(rest == 1);
        if (n <= 1000000000000LL) assert(is_prime(n) == brute_is_prime(n));
    }

    /// Strong pseudoprimes to bases 2..7 and 2..31 (151 * 751 * 28351, 149491 * 747451 * 34233211),
    /// random inputs never hit them, and only the full base set rejects them
    assert(!is_prime(3215031751LL));
    assert(!is_prime(3825123056546413051LL));

    for (long long m = 1; m <= 400; m++){
        for (long long a = -m; a < 2 * m; a++) assert(multiplicative_order(a, m) == brute_order((a % m + m) % m, m));
    }
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long m = stress::rand_int(1, 100000), a = stress::rand_int(0, m - 1);
        assert(multiplicative_order(a, m) == brute_order(a, m));
    }

    for (long long m = 1; m <= 2000; m++) assert(primitive_root(m) == brute_primitive_root(m));

    for (long long p : primes){
        vector<long long> smallest(p, -1);
        for (long long x = p - 1; x >= 0; x--) smallest[x * x % p] = x;
        for (long long a = -p; a < 2 * p; a++) assert(sqrt_mod(a, p) == smallest[(a % p + p) % p]);
    }
    for (long long p : primes){
        if (p >= 100) break;
        for (long long k = 0; k <= 2 * (p - 1); k++) check_kth_roots(k, p);
    }
    for (long long it = 0; it < stress::scaled(200); it++){
        long long p = primes[stress::rand_int(0, primes.size() - 1)];
        check_kth_roots(stress::rand_int(3, 3 * p), p);
        check_kth_roots(stress::rand_int(1, 1000000000000000000LL), p);
    }

    /// Random 62-bit primes, the largest prime below 2^62, and two with a deep 2-adic p - 1 for Tonelli-Shanks
    /// k <= 10^6 keeps the baby step giant step tables small
    const vector<long long> fixed_primes = {998244353, 29 * (1LL << 57) + 1, (1LL << 62) - 57};
    for (long long it = 0; it < stress::scaled(300); it++){
        long long p = it % 3 ? random_prime(1LL << 61, (1LL << 62) - 1) : fixed_primes[it / 3 % 3];
        long long x = stress::rand_int(1, p - 1), a = mul_mod(x, x, p), r = sqrt_mod(a, p);
        assert(r == min(x, p - x));
        long long n = stress::rand_int(1, p - 1);
        assert((sqrt_mod(n, p) == -1) == (pow_mod(n, (p - 1) / 2, p) == p - 1));

        long long order = multiplicative_order(x, p);
        assert((p - 1) % order == 0 && pow_mod(x, order, p) == 1);
        for (auto [q, e] : factorize(order)) assert(pow_mod(x, order / q, p) != 1);

        long long g = primitive_root(p);
        assert(multiplicative_order(g, p) == p - 1);
        for (long long h = 1; h < g; h++) assert(multiplicative_order(h, p) < p - 1);

        long long k = stress::rand_int(1, 1000000);
        long long y = kth_root(pow_mod(x, k, p), k, p);
        assert(pow_mod(y, k, p) == pow_mod(x, k, p));
        if (gcd(k, p - 1) <= 5000){
            vector<long long> roots = kth_roots(pow_mod(x, k, p), k, p);
            assert((long long)roots.size() == gcd(k, p - 1) && adjacent_find(roots.begin(), roots.end()) == roots.end());
            assert(binary_search(roots.begin(), roots.end(), x));
            for (long long z : roots) assert(pow_mod(z, k, p) == pow_mod(x, k, p));
        }
    }

    /// gcd(k, p - 1) holds the full power of a huge prime q in p - 1, so every discrete log is trivial:
    /// a baby step table of sqrt(q) entries would take 16 GB for the safe prime and 0.1 s per call for 2^20 q + 1
    const long long safe = 2305843009213699919LL, deep = 1048576000185597953LL;
    const long long safe_q = (safe - 1) / 2, deep_q = (deep - 1) >> 20;
    long long one_root = kth_root(1, safe - 1, safe);
    assert(1 <= one_root && one_root < safe && pow_mod(one_root, safe - 1, safe) == 1);
    for (long long it = 0; it < stress::scaled(300); it++){
        long long x = stress::rand_int(1, safe - 1), a = pow_mod(x, safe_q, safe), y = kth_root(a, safe_q, safe);
        assert(0 < y && y < safe && pow_mod(y, safe_q, safe) == a);

        long long k = deep_q * stress::rand_int(1, 1000000);
        x = stress::rand_int(1, deep - 1), a = pow_mod(x, k, deep), y = kth_root(a, k, deep);
        assert(0 < y && y < deep && pow_mod(y, k, deep) == a);
    }

    return 0;
}
