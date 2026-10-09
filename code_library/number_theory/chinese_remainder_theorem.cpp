/***
 *
 * Finds the unique solution x modulo M (product of mods) for which x % mods[i] = rems[i]
 * Mods must be pairwise co-prime
 * Product of mods must fit in int64 and each mod must be below 3 * 10^9, rems can be negative
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

int64_t mod_inverse(int64_t a, int64_t m){
    int64_t x, y, inv;
    exgcd(a, m, x, y);
    inv = (x + m) % m;
    return inv;
}

int64_t CRT(const vector<int64_t>& rems, const vector<int64_t>& mods){
    int64_t x, y, t, res = 0, prod = 1;
    for (auto mod: mods) prod *= mod;

    for (uint32_t i = 0; i < rems.size(); i++){
        x = prod / mods[i];
        y = mod_inverse(x, mods[i]);

        /// x * t mod prod only depends on t mod mods[i], so reduce first to keep every product within int64
        t = (rems[i] % mods[i] + mods[i]) % mods[i] * y % mods[i];
        res += t * x - prod;
        if (res < 0) res += prod;
    }

    return res;
}

int main(){
    auto rems = vector<int64_t>({2, 3, 2});
    auto mods = vector<int64_t>({3, 5, 7});

    assert(CRT(rems, mods) == 23);
    return 0;
}
