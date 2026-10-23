/***
 *
 * Chinese Remainder Theorem
 * Solves the system x % mods[i] = rems[i] for arbitrary positive moduli, co-prime or not
 *
 * Complexity, for k congruences and M the largest modulus:
 *   - CRT: O(k log M)
 *   - garner and the Garner constructor: O(k^2 + k log M)
 *   - g(rems): O(k^2)
 *
 * CRT(rems, mods): the unique solution x in [0, lcm(mods)), -1 if the system is inconsistent, 0 if it is empty
 *     lcm(mods) < 2^63, each mod positive, rems any long long
 * garner(rems, mods, mod): x % mod for the solution x modulo the product of the pairwise co-prime mods
 *     Each mod and the target mod in [1, 2^63), rems any long long, the product of the mods may exceed 64 bits
 * Garner g(mods, mod): garner over fixed mods with the inverses computed once, like one g(rems) per NTT coefficient
 * g(rems): garner(rems, mods, mod) with no inverse and no allocation
 *
 * Notes:
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY extended_gcd from code_library/number_theory/extended_gcd.cpp
long long extended_gcd(long long a, long long b, long long& x, long long& y){
    assert(a >= 0 && b >= 0);

    long long x1 = 0, y1 = 1;
    x = 1, y = 0;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }
    return a;
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    long long x, y;
    if (extended_gcd(a, m, x, y) != 1) return -1;
    x %= m;
    return x < 0 ? x + m : x;
}
/// END COPY extended_gcd

/// x % m + m would overflow for m above 2^62
long long normalize(long long x, long long m){
    x %= m;
    return x < 0 ? x + m : x;
}

long long CRT(const vector<long long>& rems, const vector<long long>& mods){
    long long res = 0, lcm = 1;

    for (size_t i = 0; i < mods.size(); i++){
        long long r = normalize(rems[i], mods[i]);
        long long g = gcd(lcm, mods[i]), step = mods[i] / g;
        if ((r - res) % g) return -1;


        /// res + k * lcm hits r modulo mods[i] exactly when k * (lcm / g) = (r - res) / g modulo step
        long long k = (__int128)normalize((r - res) / g, step) * mod_inverse(lcm / g % step, step) % step;
        res += k * lcm;
        lcm *= step;
    }

    return res;
}

/// Garner's digits: x = digit[0] + digit[1] * mods[0] + digit[2] * mods[0] * mods[1] + ..., m[k] is the target mod
struct Garner{
    int k;
    vector<long long> m, inv, coef, digit;

    /// coef[i * (k + 1) + j] = mods[0] * ... * mods[j - 1] % m[i], inv[i] = (mods[0] * ... * mods[i - 1])^-1 % m[i]
    Garner(const vector<long long>& mods, long long mod) : k(mods.size()), m(mods), inv(k), coef((k + 1) * (k + 1)), digit(k){
        m.push_back(mod);

        for (int i = 0; i <= k; i++){
            long long* c = &coef[i * (k + 1)];
            c[0] = 1 % m[i];
            for (int j = 1; j <= i; j++) c[j] = (unsigned __int128)c[j - 1] * m[j - 1] % m[i];
            if (i < k) inv[i] = mod_inverse(c[i], m[i]);
        }
    }

    long long operator()(const vector<long long>& rems){
        for (int i = 0; i < k; i++){
            long long diff = normalize(normalize(rems[i], m[i]) - prefix(i), m[i]);
            digit[i] = (unsigned __int128)diff * inv[i] % m[i];
        }

        return prefix(k);
    }

    /// The number written by the first i digits, % m[i]
    long long prefix(int i) const{
        const long long* c = &coef[i * (k + 1)];
        long long res = 0;
        for (int j = 0; j < i; j++) res = (res + (unsigned __int128)digit[j] * c[j]) % m[i];
        return res;
    }
};

long long garner(const vector<long long>& rems, const vector<long long>& mods, long long mod){
    return Garner(mods, mod)(rems);
}

int main(){
    const long long INF = numeric_limits<long long>::max();

    assert(CRT({2, 3, 2}, {3, 5, 7}) == 23);
    vector<long long> rems = {2, 3, 2}, mods = {3, 5, 7};
    assert(CRT(rems, mods) == 23);
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

    vector<long long> big_rems = {-1, 100000000000000000, 1LL << 62, INT64_MIN};
    vector<long long> big_mods = {1000000000000000003, 1000000000000000009, (1LL << 61) - 1, INF};
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
