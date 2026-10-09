/***
 *
 * Modular Roots
 * Multiplicative order, primitive root, square root and k-th roots modulo a prime
 *
 * Complexity:
 *   multiplicative_order(a, m): O(m^(1/4) log m) expected to factor m and phi(m), then O(log^2 m)
 *   primitive_root(m): the same factoring, then O(log^2 m) per candidate tried
 *   sqrt_mod(a, p): O(log^2 p) expected (Tonelli-Shanks)
 *   kth_root(a, k, p): O(p^(1/4) log p) expected to factor gcd(k, p - 1), then for each prime power q^e of it
 *                      O(v log p) where v is the exponent of q in p - 1 (Adleman-Manders-Miller),
 *                      plus O(v sqrt(q) log q) for baby step giant step when e < v
 *   kth_roots(a, k, p): kth_root, plus O(g log g) to list and sort the g = gcd(k, p - 1) roots
 *
 * multiplicative_order(a, m): smallest x >= 1 with a^x = 1 mod m, -1 if gcd(a, m) != 1
 * primitive_root(m): smallest g whose order is phi(m), -1 if none exists
 *                    one exists only for m = 1, 2, 4, p^k and 2p^k with p an odd prime, primitive_root(1) = 0
 * sqrt_mod(a, p): the smaller of the two x with x^2 = a mod p (the other is p - x), -1 if there is none
 * kth_root(a, k, p): any x with x^k = a mod p, -1 if there is none
 * kth_roots(a, k, p): every such x in increasing order, gcd(k, p - 1) of them when a != 0 has a root
 *
 * p must be prime (not checked), every modulus in [1, 2^62), k >= 0, a is any long long and is reduced mod m
 * 0^0 = 1, so kth_roots(1, 0, p) returns all p residues
 * kth_root needs O(sqrt(q)) memory for a prime q whose exponent in gcd(k, p - 1) is below its exponent in p - 1,
 * so a huge such q is infeasible, a huge q with equal exponents (k = q, p = 2q + 1) needs no table
 * Requires __int128 (64-bit GCC or Clang)
 *
 * The factorization is a small Pollard rho copy, not the tuned one in pollard_rho.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long mul_mod(long long a, long long b, long long m){
    return (unsigned __int128)a * b % m;
}

long long pow_mod(long long x, long long e, long long m){
    long long res = 1 % m;
    for (x %= m; e; e >>= 1, x = mul_mod(x, x, m)){
        if (e & 1) res = mul_mod(res, x, m);
    }
    return res;
}

/// Requires gcd(a, m) = 1, returns 0 for m = 1
long long inverse_mod(long long a, long long m){
    long long r = a % m, next_r = m, s = 1, next_s = 0;

    while (next_r){
        long long q = r / next_r;
        r -= q * next_r, swap(r, next_r);
        s -= q * next_s, swap(s, next_s);
    }

    return (s % m + m) % m;
}

bool is_prime(long long n){
    static const long long bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (long long q : bases){
        if (n % q == 0) return n == q;
    }

    long long d = n - 1;
    int s = 0;
    for (; !(d & 1); d >>= 1) s++;

    for (long long base : bases){
        long long x = pow_mod(base, d, n);
        if (x == 1) continue;
        for (int i = 1; i < s && x != n - 1; i++) x = mul_mod(x, x, n);
        if (x != n - 1) return false;
    }
    return true;
}

/// Brent's cycle detection with one gcd per 128 steps, returns a proper divisor of a composite n
long long pollard_rho(long long n){
    if (n % 2 == 0) return 2;

    for (long long c = 1; ; c++){
        auto next = [&](long long x){ return (mul_mod(x, x, n) + c) % n; };
        long long x = 0, y = 0, ys = 0, q = 1, g = 1;

        for (long long r = 1; g == 1; r <<= 1){
            x = y;
            for (long long i = 0; i < r; i++) y = next(y);
            for (long long k = 0; k < r && g == 1; k += 128){
                ys = y;
                for (long long i = 0; i < 128 && i < r - k; i++){
                    y = next(y);
                    q = mul_mod(q, abs(x - y), n);
                }
                g = gcd(q, n);
            }
        }

        /// The batch overshot to a multiple of n, replay it one step at a time
        if (g == n){
            do{
                ys = next(ys);
                g = gcd(abs(x - ys), n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}

void collect_primes(long long n, vector<long long>& primes){
    if (n == 1) return;
    if (is_prime(n)){
        primes.push_back(n);
        return;
    }

    long long d = pollard_rho(n);
    collect_primes(d, primes);
    collect_primes(n / d, primes);
}

/// Sorted (prime, exponent) pairs
vector<pair<long long, int>> factorize(long long n){
    vector<long long> primes;
    collect_primes(n, primes);
    sort(primes.begin(), primes.end());

    vector<pair<long long, int>> res;
    for (long long q : primes){
        if (!res.empty() && res.back().first == q) res.back().second++;
        else res.emplace_back(q, 1);
    }
    return res;
}

long long multiplicative_order(long long a, long long m){
    a = (a % m + m) % m;
    if (gcd(a, m) != 1) return -1;

    long long phi = m;
    for (auto [q, e] : factorize(m)) phi = phi / q * (q - 1);

    long long order = phi;
    for (auto [q, e] : factorize(phi)){
        while (order % q == 0 && pow_mod(a, order / q, m) == 1) order /= q;
    }
    return order;
}

long long primitive_root(long long m){
    if (m <= 4) return m - 1;

    long long odd = m % 2 ? m : m / 2;
    auto f = factorize(odd);
    if (m % 4 == 0 || f.size() != 1) return -1;

    long long phi = odd / f[0].first * (f[0].first - 1);
    auto phi_factors = factorize(phi);
    for (long long g = 2; ; g++){
        if (gcd(g, m) != 1) continue;
        bool is_root = true;
        for (auto [q, e] : phi_factors) is_root = is_root && pow_mod(g, phi / q, m) != 1;
        if (is_root) return g;
    }
}

/// Smallest n that is not a q-th power mod p, for a prime q dividing p - 1
long long non_residue(long long q, long long p){
    long long n = 2;
    while (pow_mod(n, (p - 1) / q, p) == 1) n++;
    return n;
}

long long sqrt_mod(long long a, long long p){
    a = (a % p + p) % p;
    if (a == 0 || p == 2) return a;
    if (pow_mod(a, (p - 1) / 2, p) != 1) return -1;

    long long s = p - 1;
    int r = 0;
    for (; s % 2 == 0; s /= 2) r++;

    /// Invariant: x^2 = a * b, and the order of b is 2^m for some m < r
    long long x = pow_mod(a, (s + 1) / 2, p), b = pow_mod(a, s, p), g = pow_mod(non_residue(2, p), s, p);
    while (b != 1){
        int m = 0;
        for (long long t = b; t != 1; t = mul_mod(t, t, p)) m++;

        long long gs = pow_mod(g, 1LL << (r - m - 1), p);
        g = mul_mod(gs, gs, p), x = mul_mod(x, gs, p), b = mul_mod(b, g, p), r = m;
    }

    return min(x, p - x);
}

/// Some x with x^(q^e) = b mod p, for a prime power q^e dividing p - 1 and b != 0 a q^e-th power
long long prime_power_root(long long b, long long q, int e, long long p){
    long long t = p - 1, qe = 1;
    int s = 0;
    for (; t % q == 0; t /= q) s++;
    for (int i = 0; i < e; i++) qe *= q;

    /// x^(q^e) / b = (b^t)^((q^e - 1) / t mod q^e) lies in the subgroup of order q^s, which z generates
    long long d = ((unsigned __int128)(qe - 1) * inverse_mod(t % qe, qe) % qe * t + 1) / qe;
    long long x = pow_mod(b, d, p), err = mul_mod(b, inverse_mod(pow_mod(x, qe, p), p), p);
    long long z = pow_mod(non_residue(q, p), t, p), w = pow_mod(z, (p - 1) / t / q, p);

    long long m = sqrtl(q) + 1, giant = 0;
    vector<pair<long long, long long>> baby;

    /// Baby step giant step in the subgroup of order q generated by w
    /// The table is built on the first h != 1: when q^e is the full power of q in p - 1, every h is 1
    auto log_w = [&](long long h){
        if (h == 1) return 0LL;
        if (baby.empty()){
            baby.resize(m);
            for (long long j = 0, cur = 1; j < m; j++, cur = mul_mod(cur, w, p)) baby[j] = {cur, j};
            sort(baby.begin(), baby.end());
            giant = inverse_mod(pow_mod(w, m, p), p);
        }

        for (long long i = 0; ; i++, h = mul_mod(h, giant, p)){
            auto it = lower_bound(baby.begin(), baby.end(), make_pair(h, 0LL));
            if (it != baby.end() && it->first == h) return i * m + it->second;
        }
    };

    /// Pohlig-Hellman: err = z^log_err, one base q digit per step, the lowest e digits come out 0
    long long z_inv = inverse_mod(z, p), log_err = 0;
    for (long long i = 0, qi = 1; i < s; i++, qi *= q){
        long long h = mul_mod(err, pow_mod(z_inv, log_err, p), p);
        log_err += log_w(pow_mod(h, (p - 1) / t / qi / q, p)) * qi;
    }

    return mul_mod(x, pow_mod(z, log_err / qe, p), p);
}

long long kth_root(long long a, long long k, long long p){
    a = (a % p + p) % p;
    if (k == 0) return a == 1 ? 0 : -1;
    if (a == 0) return 0;

    long long g = gcd(k, p - 1), m = (p - 1) / g;
    if (pow_mod(a, m, p) != 1) return -1;

    /// a is a g-th power, so b = a^(1 / (k / g)) is one too, and any g-th root of b is a k-th root of a
    long long b = pow_mod(a, inverse_mod(k / g % m, m), p);
    for (auto [q, e] : factorize(g)) b = prime_power_root(b, q, e, p);
    return b;
}

vector<long long> kth_roots(long long a, long long k, long long p){
    a = (a % p + p) % p;
    if (k == 0 && a == 1){
        vector<long long> all(p);
        iota(all.begin(), all.end(), 0);
        return all;
    }

    long long x = kth_root(a, k, p);
    if (x <= 0) return x ? vector<long long>{} : vector<long long>{0};

    /// w has order exactly g, so x * w^i for i < g are the g distinct roots
    long long g = gcd(k, p - 1), w = 1;
    for (auto [q, e] : factorize(g)){
        long long cofactor = p - 1;
        for (int i = 0; i < e; i++) cofactor /= q;
        w = mul_mod(w, pow_mod(non_residue(q, p), cofactor, p), p);
    }

    vector<long long> roots(g);
    for (long long i = 0; i < g; i++, x = mul_mod(x, w, p)) roots[i] = x;
    sort(roots.begin(), roots.end());
    return roots;
}

int main(){
    assert(multiplicative_order(2, 7) == 3);   /// 2^3 = 8 = 1 mod 7
    assert(multiplicative_order(4, 6) == -1);  /// gcd(4, 6) != 1
    assert(primitive_root(7) == 3);
    assert(primitive_root(998244353) == 3);
    assert(primitive_root(8) == -1);           /// none mod 2^k for k >= 3

    assert(sqrt_mod(2, 7) == 3);               /// 3^2 = 4^2 = 2 mod 7, the smaller is returned
    assert(sqrt_mod(3, 7) == -1);
    assert(sqrt_mod(-1, 13) == 5);             /// 5^2 = 25 = -1 mod 13
    assert(sqrt_mod(1000000000000000000LL, (1LL << 62) - 57) == 1000000000);

    long long x = kth_root(8, 3, 13);
    assert(pow_mod(x, 3, 13) == 8);
    assert((kth_roots(8, 3, 13) == vector<long long>{2, 5, 6}));
    assert((kth_roots(4, 4, 13) == vector<long long>{}));  /// 4 is not a 4th power mod 13
    return 0;
}
