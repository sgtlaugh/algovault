#include "common.h"

#define main library_main
#include "../code_library/power_tower.cpp"
#undef main

/// Exact value of a tower capped at CAP, by plain repeated multiplication
long long capped_tower(const vector<long long>& a, int i, long long cap){
    if (i + 1 == (int)a.size()) return min(a[i], cap);
    long long e = capped_tower(a, i + 1, cap);
    if (a[i] == 1) return 1;
    __int128 res = 1;
    for (long long k = 0; k < e; k++){
        res *= a[i];
        if (res >= cap) return cap;
    }
    return (long long)res;
}

/// Tower modulo m by cycle detection on powers of the base, no Euler theorem involved
long long brute(const vector<long long>& a, int i, long long m){
    if (m == 1) return 0;
    if (i + 1 == (int)a.size()) return a[i] % m;

    /// Powers b^0, b^1, ... mod m enter a cycle after mu steps with period lambda
    long long b = a[i] % m;
    vector<long long> first_seen(m, -1), powers;
    long long x = 1 % m, k = 0;
    while (first_seen[x] == -1){
        first_seen[x] = k++, powers.push_back(x);
        x = x * b % m;
    }
    long long mu = first_seen[x], lambda = k - mu;

    long long e = capped_tower(a, i + 1, mu + lambda + 1);
    if (e < mu + lambda) return powers[e];
    long long r = brute(a, i + 1, lambda);
    long long idx = mu + ((r - mu) % lambda + lambda) % lambda;
    return powers[idx];
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 6);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(1, it % 2 ? 4 : 12);
        long long m = stress::rand_int(1, 300);
        assert(power_tower(a, m) == brute(a, 0, m));
    }

    /// Towers whose exact value fits, checked directly
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 4);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(1, 5);
        long long exact = capped_tower(a, 0, PowerTower::LIMIT);
        if (exact >= PowerTower::LIMIT) continue;
        long long m = stress::rand_int(0, 1) ? stress::rand_int(1, 1000000000000LL) : stress::rand_int(1, 1000);
        assert(power_tower(a, m) == exact % m);
    }

    /// Range queries against whole-tower evaluation, large bases and moduli at the 1e12 bound
    for (long long it = 0; it < stress::scaled(300); it++){
        long long m = it % 3 == 0 ? 1000000000000LL : stress::rand_int(1, 1000000000000LL);
        PowerTower tower(m);
        vector<long long> a(30);
        for (auto& x : a) x = stress::rand_int(0, 3) ? stress::rand_int(1, 1000000000000000000LL) : stress::rand_int(1, 3);
        for (int q = 0; q < 20; q++){
            int l = stress::rand_int(0, 29), r = stress::rand_int(l, 29);
            vector<long long> part(a.begin() + l, a.begin() + r + 1);
            assert(tower.query(a, l, r) == power_tower(part, m));
            if (m <= 300) assert(tower.query(a, l, r) == brute(part, 0, m));
        }
    }

    vector<long long> ones(100000, 1);
    ones[0] = 7;
    assert(power_tower(ones, 1000) == 7);
    vector<long long> twos(100000, 2);
    assert(power_tower(twos, 1000000000000LL) == power_tower(vector<long long>(twos.begin(), twos.begin() + 60), 1000000000000LL));
    return 0;
}
