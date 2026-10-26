/***
 *
 * Discrete Logarithm
 * Generalized discrete logarithm with Shanks' baby step giant step, and Pohlig-Hellman when mod is prime
 *
 * Complexity:
 *   composite mod: O(sqrt(mod)) expected
 *   prime mod: O(sqrt(mod)) to factor mod - 1, then O(e sqrt(q)) expected per prime power q^e of the order of a
 *
 * discrete_log(a, b, mod): the smallest x >= 0 with a^x = b (modulo mod), -1 if no solution exists
 *
 * Notes:
 *   - a and b are any int (negative ones too), 1 <= mod <= INT_MAX, a, b and mod need not be pairwise co-prime
 *   - a^0 = 1 (including 0^0), so x = 0 whenever b = 1 (modulo mod), which is every b when mod = 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// The table masks the low bits, and residues like b * 3^i mod 2^30 all share them
struct mix_hash{
    size_t operator()(int x) const{
        return (unsigned long long)x * 0x9E3779B97F4A7C15ULL >> 32;
    }
};

namespace dlog{
    int pow_mod(long long x, long long k, int mod){
        long long res = 1 % mod;
        for (; k; k >>= 1, x = x * x % mod){
            if (k & 1) res = res * x % mod;
        }
        return res;
    }

    /// Miller-Rabin with the bases 2, 7 and 61, deterministic below 4759123141 (Jaeschke 1993)
    bool is_prime(int n){
        if (n < 2 || n % 2 == 0) return n == 2;

        int s = __builtin_ctz(n - 1);
        for (int w : {2, 7, 61}){
            if (w % n == 0) return true;
            long long x = pow_mod(w, (n - 1) >> s, n);
            if (x == 1) continue;
            for (int r = 1; r < s && x != n - 1; r++) x = x * x % n;
            if (x != n - 1) return false;
        }

        return true;
    }

    /// Smallest x >= 0 with a^x = b for 0 <= a, b < mod, given that any solution has one below bound
    int baby_step_giant_step(int a, int b, int mod, int bound){
        int n = sqrt(bound + 0.5) + 1;

        /// Linear probing at most half full, -1 marks empty: 51% faster than gp_hash_table on Library Checker data
        int cap = 2 << __lg(2 * n + 5);
        vector<int> keys(cap, -1), vals(cap);
        auto slot = [&](int key){
            size_t h = mix_hash()(key) & (cap - 1);
            while (keys[h] != -1 && keys[h] != key) h = (h + 1) & (cap - 1);
            return h;
        };

        int cur = 1 % mod;
        for (int i = 0; i < n + 3; i++){
            if (cur == b) return i;
            int key = (long long)b * cur % mod;
            size_t h = slot(key);
            keys[h] = key, vals[h] = i;
            cur = (long long)cur * a % mod;
        }

        int giant = pow_mod(a, n, mod);
        cur = giant;
        for (int i = 2; i < n + 3; i++){
            cur = (long long)cur * giant % mod;
            size_t h = slot(cur);
            if (keys[h] == cur){
                int x = (((long long)n * i - vals[h]) % mod + mod) % mod;
                return pow_mod(a, x, mod) == b ? x : -1;
            }
        }

        return -1;
    }

    /// Smallest x >= 0 with a^x = b for a prime p and 0 < a < p, 0 <= b < p
    int pohlig_hellman(int a, int b, int p){
        vector<int> primes;
        int rest = p - 1;
        for (int q = 2; q <= rest / q; q++){
            if (rest % q) continue;
            primes.push_back(q);
            while (rest % q == 0) rest /= q;
        }
        if (rest > 1) primes.push_back(rest);

        int ord = p - 1;
        for (int q : primes){
            while (ord % q == 0 && pow_mod(a, ord / q, p) == 1) ord /= q;
        }
        /// (Z/p)* is cyclic, so the powers of a are exactly the y with y^ord = 1
        if (pow_mod(b, ord, p) != 1) return -1;

        /// x is the answer modulo m, the product of the prime powers of ord handled so far
        long long x = 0, m = 1;
        for (int q : primes){
            long long qe = 1;
            while (ord % (qe * q) == 0) qe *= q;
            if (qe == 1) continue;

            /// g has order q^e, so x mod q^e is found one base q digit at a time with gamma of order q
            int g = pow_mod(a, ord / qe, p), h = pow_mod(b, ord / qe, p);
            int gamma = pow_mod(g, qe / q, p), g_inv = pow_mod(g, qe - 1, p);
            long long x_qe = 0;
            for (long long qk = 1; qk < qe; qk *= q){
                int hk = pow_mod((long long)h * pow_mod(g_inv, x_qe, p) % p, qe / qk / q, p);
                x_qe += baby_step_giant_step(gamma, hk, p, q) * qk;
            }

            /// Euler's theorem: m^-1 = m^(phi(q^e) - 1) modulo q^e
            long long m_inv = pow_mod(m % qe, qe / q * (q - 1) - 1, qe);
            x += m * ((x_qe - x % qe + qe) % qe * m_inv % qe);
            m *= qe;
        }

        return x;
    }
}

int discrete_log(int a, int b, int mod){
    a %= mod, b %= mod;
    if (a < 0) a += mod;
    if (b < 0) b += mod;

    if (a && dlog::is_prime(mod)) return dlog::pohlig_hellman(a, b, mod);
    return dlog::baby_step_giant_step(a, b, mod, mod);
}

int main(){
    /// Smallest x with a^x = b, a, b and mod need not be co-prime
    assert(discrete_log(2, 6, 10) == 4);       /// 2^4 = 16 = 6 mod 10
    assert(discrete_log(5, 33, 58) == 9);
    assert(discrete_log(2, 3, 4) == -1);       /// powers of 2 mod 4 are 1, 2, 0, 0, ...
    assert(discrete_log(2, 1, 3) == 0);        /// a^0 = 1
    assert(discrete_log(-2, 3, 5) == 1);       /// -2 = 3 mod 5
    return 0;
}
