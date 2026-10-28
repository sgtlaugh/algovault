#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/discrete_logarithm.cpp"
#undef main

int residue(long long v, int mod){
    return (v % mod + mod) % mod;
}

int pow_mod(long long x, int n, int mod){
    long long res = 1 % mod;
    for (x = residue(x, mod); n; n >>= 1, x = x * x % mod){
        if (n & 1) res = res * x % mod;
    }

    return res;
}

int brute(int a, int b, int mod){
    a = residue(a, mod), b = residue(b, mod);
    long long e = 1 % mod;
    for (int x = 0; x <= mod; x++, e = e * a % mod){
        if (e == b) return x;
    }

    return -1;
}

bool slow_is_prime(long long n){
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; d++){
        if (n % d == 0) return false;
    }

    return true;
}

/// A random prime p <= hi, above hi / max_factor, whose p - 1 has only prime factors up to max_factor
int smooth_prime(long long hi, int max_factor){
    vector<int> small;
    for (int q = 2; q <= max_factor; q++){
        if (slow_is_prime(q)) small.push_back(q);
    }

    while (true){
        long long p = 2;
        while (p * max_factor <= hi) p *= small[stress::rand_int(0, (int)small.size() - 1)];
        if (p + 1 <= hi && slow_is_prime(p + 1)) return p + 1;
    }
}

int main(){
    /// Strong pseudoprimes to two of the bases 2, 7 and 61 are composite, n < 5000 holds those of one base
    for (int n : {79381, 314821, 916327}) assert(!dlog::is_prime(n));
    for (int n = 0; n < 5000; n++) assert(dlog::is_prime(n) == slow_is_prime(n));
    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = stress::rand_int(1, INT_MAX);
        assert(dlog::is_prime(n) == slow_is_prime(n));
    }

    /// Every a and b for small primes: the Pohlig-Hellman path, with p = 2 and a = 0 among them
    for (int p = 2; p < 150; p++){
        if (!slow_is_prime(p)) continue;
        for (int a = 0; a < p; a++){
            for (int b = 0; b < p; b++) assert(discrete_log(a, b, p) == brute(a, b, p));
        }
    }

    /// Large primes against the baby step giant step path: random, smooth p - 1, and p = 2q + 1 with q prime
    for (long long it = 0; it < stress::scaled(150); it++){
        long long p = 0, hi = it % 4 ? INT_MAX : 1000000;
        if (it % 3 == 0){
            do p = stress::rand_int(2, hi); while (!slow_is_prime(p));
        }
        else if (it % 3 == 1) p = smooth_prime(hi, stress::rand_int(0, 1) ? 50 : 3000);
        else{
            do p = 2 * stress::rand_int(1, (hi - 1) / 2) + 1; while (!slow_is_prime(p) || !slow_is_prime(p / 2));
        }

        int a = stress::rand_int(it % 5 ? 1 : INT_MIN, INT_MAX), k = stress::rand_int(0, p - 2);
        int b = stress::rand_int(0, 1) ? pow_mod(a, k, p) : stress::rand_int(INT_MIN, INT_MAX);
        assert(discrete_log(a, b, p) == dlog::baby_step_giant_step(residue(a, p), residue(b, p), p, p));
    }

    /// The largest primes, 2^31 - 1 has the prime power 3^2 in p - 1
    for (int p : {INT_MAX, 2147483629, 2147483587}){
        for (int it = 0; it < 4; it++){
            int a = stress::rand_int(1, p - 1);
            int b = it % 2 ? pow_mod(a, stress::rand_int(0, p - 2), p) : stress::rand_int(0, p - 1);
            assert(discrete_log(a, b, p) == dlog::baby_step_giant_step(a, b, p, p));
        }
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        int mod = stress::rand_int(1, it % 10 ? 100 : 5000);
        int a = stress::rand_int(it % 5 ? 0 : -3 * mod, it % 3 ? mod - 1 : 3 * mod);
        int b = it % 2 ? pow_mod(a, stress::rand_int(0, 2 * mod), mod) : stress::rand_int(-3 * mod, 3 * mod);
        if (it % 4 == 1) b -= mod * stress::rand_int(0, 3);
        assert(discrete_log(a, b, mod) == brute(a, b, mod));
    }

    /// Moduli too large for brute force: b = a^k, so a solution no larger than k must come back
    for (long long it = 0; it < stress::scaled(30); it++){
        int mod = stress::rand_int(it % 2 ? 1000000 : 2000000000, INT_MAX), a = stress::rand_int(it % 3 ? 0 : INT_MIN, mod - 1), k = stress::rand_int(0, mod);
        int b = pow_mod(a, k, mod), x = discrete_log(a, b, mod);
        assert(0 <= x && x <= k && pow_mod(a, x, mod) == b);
    }

    /// INT_MIN for both a and b, whose reduction would overflow if negated first
    assert(discrete_log(INT_MIN, INT_MIN, INT_MAX) == 1);

    /// Baby steps b * 3^i mod 2^30 share their low bits, the worst case for an identity hash with masked buckets
    auto start = chrono::steady_clock::now();
    for (int shift = 11; shift <= 14; shift++){
        for (int odd = 5; odd <= 13; odd += 2) assert(discrete_log(3, odd << shift, 1 << 30) == -1);
    }
    assert(chrono::steady_clock::now() - start < chrono::seconds(10));

    return 0;
}
