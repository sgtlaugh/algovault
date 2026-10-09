/***
 *
 * Chinese Remainder Theorem
 * Solves the system x % mods[i] = rems[i] for arbitrary positive moduli, co-prime or not
 *
 * Complexity: O(k log M) for CRT, O(k^2 + k log M) for garner and the Garner build, O(k^2) per Garner query
 *   k congruences, M the largest modulus
 * Requires __int128
 *
 * CRT(rems, mods) returns the unique solution x in [0, lcm(mods)), or -1 if the system is inconsistent
 *   The lcm of the mods must be below 2^63, each mod must be positive, rems can be any int64 (negative too)
 *   An empty system returns 0
 *
 * garner(rems, mods, mod) returns x % mod for the solution x modulo the product of the mods (Garner's mixed radix)
 *   For pairwise co-prime mods whose product does not fit in 64 bits, e.g. recombining NTT results
 *   Each mod and the target mod must be in [1, 2^63), rems can be any int64
 *
 * Garner g(mods, mod) precomputes the inverses once, then g(rems) = garner(rems, mods, mod) with no inverse and no allocation
 *   Use it when the mods repeat across calls, e.g. one g(rems) per coefficient of an NTT convolution
 *
 * Example:
 *   CRT({2, 4}, {4, 6}) = 10, CRT({1, 2}, {4, 6}) = -1
 *   garner({1, 2, 3}, {1000000007, 998244353, 1000000009}, 1000000007) = 1
 *   Garner g({998244353, 167772161, 469762049}, 1000000007); g({1, 2, 3}) = 330560083
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// Bezout's identity, ax + by = gcd(a,b)
int64_t exgcd(int64_t a, int64_t b, int64_t& x, int64_t& y){
    if (!b){
        y = 0, x = 1;
        return a;
    }

    int64_t g = exgcd(b, a % b, y, x);
    y -= (a / b) * x;
    return g;
}

/// x % m + m would overflow for m above 2^62
int64_t normalize(int64_t x, int64_t m){
    x %= m;
    return x < 0 ? x + m : x;
}

int64_t mod_inverse(int64_t a, int64_t m){
    int64_t x, y;
    exgcd(a, m, x, y);
    return normalize(x, m);
}

int64_t CRT(const vector<int64_t>& rems, const vector<int64_t>& mods){
    int64_t res = 0, lcm = 1;

    for (size_t i = 0; i < mods.size(); i++){
        int64_t r = normalize(rems[i], mods[i]);
        int64_t g = gcd(lcm, mods[i]), step = mods[i] / g;
        if ((r - res) % g) return -1;


        /// res + k * lcm hits r modulo mods[i] exactly when k * (lcm / g) = (r - res) / g modulo step
        int64_t k = (__int128)normalize((r - res) / g, step) * mod_inverse(lcm / g % step, step) % step;
        res += k * lcm;
        lcm *= step;
    }

    return res;
}

/// Garner's digits: x = digit[0] + digit[1] * mods[0] + digit[2] * mods[0] * mods[1] + ..., m[k] is the target mod
struct Garner{
    int k;
    vector<int64_t> m, inv, coef, digit;

    /// coef[i * (k + 1) + j] = mods[0] * ... * mods[j - 1] % m[i], inv[i] = (mods[0] * ... * mods[i - 1])^-1 % m[i]
    Garner(const vector<int64_t>& mods, int64_t mod) : k(mods.size()), m(mods), inv(k), coef((k + 1) * (k + 1)), digit(k){
        m.push_back(mod);

        for (int i = 0; i <= k; i++){
            int64_t* c = &coef[i * (k + 1)];
            c[0] = 1 % m[i];
            for (int j = 1; j <= i; j++) c[j] = (unsigned __int128)c[j - 1] * m[j - 1] % m[i];
            if (i < k) inv[i] = mod_inverse(c[i], m[i]);
        }
    }

    int64_t operator()(const vector<int64_t>& rems){
        for (int i = 0; i < k; i++){
            int64_t diff = normalize(normalize(rems[i], m[i]) - prefix(i), m[i]);
            digit[i] = (unsigned __int128)diff * inv[i] % m[i];
        }

        return prefix(k);
    }

    /// The number written by the first i digits, % m[i]
    int64_t prefix(int i) const{
        const int64_t* c = &coef[i * (k + 1)];
        int64_t res = 0;
        for (int j = 0; j < i; j++) res = (res + (unsigned __int128)digit[j] * c[j]) % m[i];
        return res;
    }
};

int64_t garner(const vector<int64_t>& rems, const vector<int64_t>& mods, int64_t mod){
    return Garner(mods, mod)(rems);
}

int main(){
    const int64_t INF = numeric_limits<int64_t>::max();

    assert(CRT({2, 3, 2}, {3, 5, 7}) == 23);
    assert(CRT({}, {}) == 0);
    assert(CRT({5}, {1}) == 0);
    assert(CRT({-7}, {5}) == 3);
    assert(CRT({2, 4}, {4, 6}) == 10);
    assert(CRT({3, 5}, {10, 12}) == 53);
    assert(CRT({-1, -1}, {4, 6}) == 11);
    assert(CRT({1, 2}, {4, 6}) == -1);
    assert(CRT({7, 2}, {5, 5}) == 2);
    assert(CRT({1, 2}, {5, 5}) == -1);
    assert(CRT({1, 2, 3}, {2, 3, 4}) == 11);
    assert(CRT({0, 2, 3}, {2, 3, 4}) == -1);
    assert(CRT({5, 3 * (1LL << 60) + 5}, {3 * (1LL << 61), 3 * (1LL << 60)}) == 5);
    assert(CRT({5, 6}, {3 * (1LL << 61), 3 * (1LL << 60)}) == -1);
    assert(CRT({-51146671749066063, 9680838726321}, {1099511627776LL * 59049, 1073741824LL * 1701}) == 78703452468023985);
    assert(CRT({-1}, {INF}) == INF - 1);
    assert(CRT({INT64_MIN}, {INF}) == INF - 1);
    assert(CRT({0, 1}, {153092023, 60247241209}) == 2695413193155969174);
    assert(CRT({INT64_MIN, INF}, {153092023, 60247241209}) == 2695413193155969173);

    assert(garner({2, 3, 2}, {3, 5, 7}, 10) == 3);
    assert(garner({}, {}, 10) == 0);
    assert(garner({2, 3, 2}, {3, 5, 7}, 1) == 0);
    assert(garner({1, 2, 3}, {1000000007, 998244353, 1000000009}, 1000000007) == 1);
    assert(garner({1, 2, 3}, {1000000007, 998244353, 1000000009}, (1LL << 61) - 1) == 1131585747248034387);
    assert(garner({1, 2, 3}, {1000000007, 998244353, 1000000009}, INF) == 1131585747206765391);

    vector<int64_t> big_rems = {-1, 100000000000000000, 1LL << 62, INT64_MIN};
    vector<int64_t> big_mods = {1000000000000000003, 1000000000000000009, (1LL << 61) - 1, INF};
    assert(garner(big_rems, big_mods, INF) == INF - 1);
    assert(garner(big_rems, big_mods, 1000000007) == 746669198);
    assert(garner(big_rems, big_mods, 2) == 0);

    Garner ntt({998244353, 167772161, 469762049}, 1000000007);
    assert(ntt({1, 2, 3}) == 330560083);
    assert(ntt({0, 0, 0}) == 0);
    assert(ntt({1, 2, 3}) == 330560083);

    Garner empty({}, 10);
    assert(empty({}) == 0);

    return 0;
}
