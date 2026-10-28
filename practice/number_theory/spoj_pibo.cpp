/***
 *
 * Chinese Remainder Theorem
 * Solves the system x % mods[i] = rems[i] for arbitrary positive moduli, co-prime or not
 *
 * Complexity, for k congruences and M the largest modulus:
 *   - CRT: O(k log M)
 *   - Garner g(mods, mod): O(k^2 + k log M)
 *   - g(rems): O(k^2)
 *
 * CRT(rems, mods): the unique solution x in [0, lcm(mods)), -1 if the system is inconsistent, 0 if it is empty
 * Garner g(mods, mod): prepares g(rems) for pairwise co-prime mods, asserted, computing the inverses once
 * g(rems): x % mod for the solution x in [0, product of the mods), allocation free, like one call per NTT coefficient
 *
 * Notes:
 *   - CRT: each mod positive, lcm(mods) < 2^63 asserted, rems any long long
 *   - Garner: each mod and the target mod in [1, 2^63), rems any long long, the product of the mods may exceed 64 bits
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY extended_gcd from code_library/number_theory/extended_gcd.cpp
array<long long, 3> extended_gcd(long long a, long long b){
    assert(a >= 0 && b >= 0);

    long long x = 1, y = 0, x1 = 0, y1 = 1;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }

    return {a, x, y};
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    auto [g, x, y] = extended_gcd(a, m);
    if (g != 1) return -1;
    return x < 0 ? x + m : x;
}
/// END COPY extended_gcd

/// x % m + m would overflow for m above 2^62
long long non_negative_mod(long long x, long long m){
    x %= m;
    return x < 0 ? x + m : x;
}

long long CRT(const vector<long long>& rems, const vector<long long>& mods){
    long long res = 0, lcm = 1;

    for (size_t i = 0; i < mods.size(); i++){
        long long r = non_negative_mod(rems[i], mods[i]);
        long long g = gcd(lcm, mods[i]), step = mods[i] / g;
        if ((r - res) % g) return -1;
        assert((__int128)lcm * step <= LLONG_MAX);

        /// res + k * lcm hits r modulo mods[i] exactly when k * (lcm / g) = (r - res) / g modulo step
        long long k = (__int128)non_negative_mod((r - res) / g, step) * mod_inverse(lcm / g, step) % step;
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
            c[0] = 1;
            for (int j = 1; j <= i; j++) c[j] = (unsigned __int128)c[j - 1] * m[j - 1] % m[i];
            if (i < k){
                inv[i] = mod_inverse(c[i], m[i]);
                assert(inv[i] != -1);
            }
        }
    }

    long long operator()(const vector<long long>& rems){
        for (int i = 0; i < k; i++){
            long long diff = non_negative_mod(rems[i], m[i]) - prefix(i);
            if (diff < 0) diff += m[i];
            digit[i] = (unsigned __int128)diff * inv[i] % m[i];
        }

        return prefix(k);
    }

    /// The number written by the first i digits, % m[i]
    long long prefix(int i) const{
        const long long* c = &coef[i * (k + 1)];
        unsigned __int128 res = 0;
        for (int j = 0; j < i; j++){
            res += (unsigned __int128)digit[j] * c[j];
            if (j % 3 == 2) res %= m[i];  /// 3 products below 2^126 plus a residue fit 128 bits, 7-36% faster
        }

        return res % m[i];
    }
};

#define MAX 105
#define MAXP 9100

int n, d, c[MAX], val[MAXP], pre[MAXP], fib[MAXP];

long long range_sum(int p, int s, long long len){
    int r = len % p;
    long long res = (len / p) % p * pre[p];

    if ((s + r) <= p) res += pre[s + r] - pre[s];
    else res += pre[p] - pre[s] + pre[s + r - p];
    return res % p;
}

long long solve(int p){
    int i, j;
    long long res;

    for (i = 0; i < p; i++){
        val[i] = 0;
        for (j = d; j >= 0; j--) val[i] = (val[i] * i + c[j]) % p;
        pre[i + 1] = (pre[i] + val[i]) % p;
    }

    fib[0] = 0, fib[1] = 1;
    for (i = 2; i < p; i++) fib[i] = (fib[i - 1] + fib[i - 2]) % p;

    res = fib[(n + 1) % (p - 1)];
    for (i = 1; i < p && i < n; i++){
        res += fib[i] * range_sum(p, (n + 1 - i) % p, (n - 1 - i) / (p - 1) + 1);
    }

    return res % p;
}

int main(){
    int i;

    scanf("%d %d", &n, &d);
    for (i = 0; i <= d; i++) scanf("%d", &c[i]);

    printf("%lld\n", CRT({solve(11), solve(41), solve(271), solve(9091)}, {11, 41, 271, 9091}));
    return 0;
}
