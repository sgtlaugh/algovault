/***
 *
 * Power Tower
 * a[l] ^ (a[l + 1] ^ (... ^ a[r])) modulo m, evaluated from the top down
 *
 * Complexity: O(sqrt(m) log m) to build, O(log^2 m) per query
 *
 * PowerTower tower(m): 1 <= m <= 1e12, precomputes m, phi(m), phi(phi(m)), ... down to 1
 * tower.query(a, l, r): the tower over a[l..r] modulo m, 1 <= a[i] <= 1e18
 * power_tower(a, m): the whole array as one tower
 *
 * Uses a^e = a^(e mod phi(m) + phi(m)) (mod m), valid once e >= phi(m) because phi(m) >= every prime exponent of m
 * Exponents below 2^62 are computed exactly instead, so the rule is only applied where it holds
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PowerTower{
    static constexpr long long LIMIT = 1LL << 62;
    vector<long long> chain;

    static long long phi(long long n){
        long long res = n;
        for (long long p = 2; p * p <= n; p++){
            if (n % p) continue;
            while (n % p == 0) n /= p;
            res = res / p * (p - 1);
        }
        return n > 1 ? res / n * (n - 1) : res;
    }

    PowerTower(long long m) : chain{m}{
        assert(1 <= m && m <= 1000000000000LL);
        while (chain.back() > 1) chain.push_back(phi(chain.back()));
    }

    static long long pow_mod(long long b, long long e, long long m){
        __int128 res = 1 % m, x = b % m;
        for (; e; e >>= 1, x = x * x % m){
            if (e & 1) res = res * x % m;
        }
        return (long long)res;
    }

    /// b^e capped at LIMIT
    static long long pow_saturated(long long b, long long e){
        if (b == 1 || e == 0) return 1;
        if (e >= 62) return LIMIT;
        __int128 res = 1;
        for (long long i = 0; i < e; i++){
            res *= b;
            if (res >= LIMIT) return LIMIT;
        }
        return (long long)res;
    }

    /// Exact value of a[i..r] capped at LIMIT, a 1 ends the tower and five values >= 2 always exceed it
    long long exact(const vector<long long>& a, int i, int r) const{
        int j = i;
        while (j < r && a[j] != 1 && j - i < 4) j++;
        long long value = min(a[j], LIMIT);
        for (int k = j - 1; k >= i; k--) value = pow_saturated(a[k], value);
        return value;
    }

    long long solve(const vector<long long>& a, int i, int r, int depth) const{
        long long m = chain[depth];
        if (m == 1) return 0;
        if (i == r || a[i] == 1) return a[i] % m;

        long long e = exact(a, i + 1, r);
        if (e < LIMIT) return pow_mod(a[i], e, m);
        long long ph = chain[depth + 1];
        return pow_mod(a[i], solve(a, i + 1, r, depth + 1) + ph, m);
    }

    long long query(const vector<long long>& a, int l, int r) const{
        assert(0 <= l && l <= r && r < (int)a.size());
        for (int i = l; i <= r; i++) assert(a[i] >= 1);
        return solve(a, l, r, 0);
    }
};

long long power_tower(const vector<long long>& a, long long m){
    return PowerTower(m).query(a, 0, a.size() - 1);
}

int main(){
    assert(power_tower({2, 3, 2}, 1000) == 512);
    assert(power_tower({3, 3, 3}, 100) == 87);
    assert(power_tower({5}, 3) == 2);
    assert(power_tower({7, 1, 100}, 10) == 7);
    assert(power_tower({123456789, 987654321}, 1) == 0);
    assert(power_tower({2, 2, 2, 2}, 100000) == 65536);
    assert(power_tower({2, 10}, 6) == 4);
    assert(power_tower({3, 2, 2}, 1000) == 81);

    PowerTower tower(1000000);
    vector<long long> a = {1, 2, 3, 4, 2};
    assert(tower.query(a, 1, 2) == 8);
    assert(tower.query(a, 0, 4) == 1);
    assert(tower.query(a, 3, 4) == 16);
    return 0;
}
