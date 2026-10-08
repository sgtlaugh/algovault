#include "common.h"

#define main library_main
#include "../code_library/pisano_period.cpp"
#undef main

typedef unsigned long long ull;

ull mul_mod(ull a, ull b, ull m){
    return (unsigned __int128)a * b % m;
}

/// Independent of the library's Montgomery arithmetic, the first 12 prime bases are exact below 3.18 * 10^23
bool reference_is_prime(ull n){
    if (n < 2) return false;
    for (ull p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) if (n % p == 0) return n == p;
    ull d = n - 1;
    int s = 0;
    while (!(d & 1)) d >>= 1, s++;
    for (ull a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        ull x = 1, b = a, e = d;
        for (; e; e >>= 1, b = mul_mod(b, b, n)) if (e & 1) x = mul_mod(x, b, n);
        bool composite = x != 1 && x != n - 1;
        for (int r = 1; r < s && composite; r++) composite = (x = mul_mod(x, x, n)) != n - 1;
        if (composite) return false;
    }
    return true;
}

/// Independent (F(k) mod m, F(k + 1) mod m) by 2x2 matrix exponentiation
pair<ull, ull> fib_pair(ull k, ull m){
    ull a = 1 % m, b = 0, c = 0, d = 1 % m;      /// result matrix, starts as the identity
    ull p = 0, q = 1 % m, r = 1 % m, s = 1 % m;  /// [[0, 1], [1, 1]]
    for (; k; k >>= 1){
        if (k & 1){
            ull na = (mul_mod(a, p, m) + mul_mod(b, r, m)) % m, nb = (mul_mod(a, q, m) + mul_mod(b, s, m)) % m;
            ull nc = (mul_mod(c, p, m) + mul_mod(d, r, m)) % m, nd = (mul_mod(c, q, m) + mul_mod(d, s, m)) % m;
            a = na, b = nb, c = nc, d = nd;
        }
        ull np = (mul_mod(p, p, m) + mul_mod(q, r, m)) % m, nq = (mul_mod(p, q, m) + mul_mod(q, s, m)) % m;
        ull nr = (mul_mod(r, p, m) + mul_mod(s, r, m)) % m, ns = (mul_mod(r, q, m) + mul_mod(s, s, m)) % m;
        p = np, q = nq, r = nr, s = ns;
    }
    return {b, d};  /// the matrix power is [[F(k - 1), F(k)], [F(k), F(k + 1)]]
}

bool is_period(ull k, ull m){
    return fib_pair(k, m) == make_pair(0 % m, 1 % m);
}

long long brute_pisano(long long n){
    if (n == 1) return 1;
    long long a = 0, b = 1, k = 0;
    do { long long c = (a + b) % n; a = b, b = c, k++; } while (!(a == 0 && b == 1));
    return k;
}

int main(){
    /// Strong pseudoprimes to the smallest bases, a weakened base set would call them prime
    for (ull n : {3215031751ULL, 2152302898747ULL, 3474749660383ULL, 341550071728321ULL, 3825123056546413051ULL}) assert(!rho::miller_rabin(n));

    rho::init();

    for (long long n = 1; n <= 3000; n++) assert(pisano_period(n) == brute_pisano(n));

    /// 299210837 is a prime that divides a Miller-Rabin base, it used to hang the factorization
    assert(rho::miller_rabin(299210837) && pisano_period(299210837) == 199473892 && pisano_period(2 * 299210837LL) == 598421676);

    /// Larger n: the answer must be a period, and no p / q for a prime q dividing it may be one
    for (long long it = 0; it < stress::scaled(60); it++){
        long long n = stress::rand_int(2, it % 2 ? 1000000000000LL : 1000000);
        long long p = pisano_period(n);
        assert(p > 0 && p <= 6 * n && is_period(p, n));
        long long rest = p;
        for (long long q = 2; q * q <= rest; q++){
            if (rest % q) continue;
            assert(!is_period(p / q, n));
            while (rest % q == 0) rest /= q;
        }
        if (rest > 1) assert(!is_period(p / rest, n));
    }

    /// The factorization and primality helpers behind it
    for (long long it = 0; it < stress::scaled(15000); it++){
        long long n = (long long)(stress::rng()() >> stress::rand_int(4, 63)) + 1;
        if (n > 1500000000000000000LL) continue;
        auto factors = rho::factorize(n);
        __int128 product = 1;
        for (auto f : factors) product *= f, assert(reference_is_prime(f));
        assert(product == n && is_sorted(factors.begin(), factors.end()));

        ull m = stress::rand_int(1, 1000000000000000LL) | 1, k = stress::rand_int(0, 1000000000);
        auto expected = fib_pair(k, m);
        assert(fib(k, m) == make_pair((long long)expected.first, (long long)expected.second));
    }
    return 0;
}
