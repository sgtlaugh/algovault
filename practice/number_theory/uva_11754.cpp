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
 *   The lcm of the mods must be below 2^63, each mod must be positive, rems can be any long long (negative too)
 *   An empty system returns 0
 *
 * garner(rems, mods, mod) returns x % mod for the solution x modulo the product of the mods (Garner's mixed radix)
 *   For pairwise co-prime mods whose product does not fit in 64 bits, e.g. recombining NTT results
 *   Each mod and the target mod must be in [1, 2^63), rems can be any long long
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
        long long g = __gcd(lcm, mods[i]), step = mods[i] / g;
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
    const long long LIMIT = 10000;
    int c, s;

    while (scanf("%d %d", &c, &s) == 2 && (c || s)){
        vector<long long> x(c);
        vector<vector<long long>> ys(c);
        long long combos = 1, product = 1;
        int best = 0;
        for (int i = 0; i < c; i++){
            int k;
            scanf("%lld %d", &x[i], &k);
            ys[i].resize(k);
            for (auto& y : ys[i]) scanf("%lld", &y);
            sort(ys[i].begin(), ys[i].end());

            product *= x[i];
            combos = min(combos * k, LIMIT + 1);
            if (ys[i].size() * x[best] < ys[best].size() * x[i]) best = i;
        }

        if (combos <= LIMIT){
            vector<long long> base, rems(c);
            function<void(int)> pick = [&](int i) -> void{
                if (i == c){
                    base.push_back(CRT(rems, x));
                    return;
                }
                for (long long y : ys[i]){
                    rems[i] = y;
                    pick(i + 1);
                }
            };
            pick(0);
            sort(base.begin(), base.end());

            for (long long d = 0; s; d++){
                for (size_t i = 0; i < base.size() && s; i++){
                    long long v = d * product + base[i];
                    if (v > 0) printf("%lld\n", v), s--;
                }
            }
        }
        else{
            vector<set<long long>> allowed(c);
            for (int i = 0; i < c; i++) allowed[i] = set<long long>(ys[i].begin(), ys[i].end());
            for (long long d = 0; s; d++){
                for (size_t i = 0; i < ys[best].size() && s; i++){
                    long long v = d * x[best] + ys[best][i];
                    bool ok = v > 0;
                    for (int j = 0; j < c && ok; j++) ok = allowed[j].count(v % x[j]);
                    if (ok) printf("%lld\n", v), s--;
                }
            }
        }
        puts("");
    }

    return 0;
}
