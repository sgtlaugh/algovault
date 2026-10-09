#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/number_theory/power_tower.cpp"
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

bool aborts(const vector<long long>& a, long long m){
    pid_t pid = fork();
    assert(pid >= 0);

    if (pid == 0){
        assert(freopen("/dev/null", "w", stderr));
        power_tower(a, m);
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    return WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}

int main(){
    /// Building a tower near 1e12 takes about a million steps, so large moduli share a prebuilt pool
    vector<PowerTower> pool;
    for (long long m : {1000000000000LL, 999999999989LL, 1LL << 39, 847288609443LL, 963761198400LL}) pool.emplace_back(m);
    for (int i = 0; i < 3; i++) pool.emplace_back(stress::rand_int(1, 1000000000000LL));

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

        if (stress::rand_int(0, 1)){
            const PowerTower& tower = pool[stress::rand_int(0, pool.size() - 1)];
            assert(tower.query(a, 0, n - 1) == exact % tower.chain[0]);
        }
        else{
            long long m = stress::rand_int(1, 1000);
            assert(power_tower(a, m) == exact % m);
        }
    }

    /// Range queries with large bases, towers far above 2^62, against the cycle-detection brute force
    vector<PowerTower> small_pool;
    for (long long m : {200000LL, 199999LL, 1LL << 17, 177147LL, 166320LL}) small_pool.emplace_back(m);
    for (int i = 0; i < 3; i++) small_pool.emplace_back(stress::rand_int(1, 200000));
    for (long long it = 0; it < stress::scaled(60); it++){
        const PowerTower& tower = small_pool[it % small_pool.size()];
        vector<long long> a(30);
        for (auto& x : a) x = stress::rand_int(0, 3) ? stress::rand_int(1, 1000000000000000000LL) : stress::rand_int(1, 3);
        for (int q = 0; q < 10; q++){
            int l = stress::rand_int(0, 29), r = stress::rand_int(l, 29);
            vector<long long> part(a.begin() + l, a.begin() + r + 1);
            assert(tower.query(a, l, r) == brute(part, 0, tower.chain[0]));
        }
    }

    /// Range queries at the 1e12 bound against the same tower on the extracted range
    for (long long it = 0; it < stress::scaled(300); it++){
        const PowerTower& tower = pool[it % pool.size()];
        vector<long long> a(30);
        for (auto& x : a) x = stress::rand_int(0, 3) ? stress::rand_int(1, 1000000000000000000LL) : stress::rand_int(1, 3);
        for (int q = 0; q < 20; q++){
            int l = stress::rand_int(0, 29), r = stress::rand_int(l, 29);
            vector<long long> part(a.begin() + l, a.begin() + r + 1);
            assert(tower.query(a, l, r) == tower.query(part, 0, r - l));
        }
    }

    /// Larger moduli, where phi chains are longer and the exact/Euler switch happens deeper
    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(1, 5);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(1, 6);
        long long m = stress::rand_int(100000, 1000000);
        assert(power_tower(a, m) == brute(a, 0, m));
    }
    assert(power_tower({2, 2, 2, 2, 2, 2}, 1009LL * 1013 * 3) == 259264);

    vector<long long> ones(100000, 1);
    ones[0] = 7;
    assert(power_tower(ones, 1000) == 7);

    vector<long long> twos(100000, 2);
    assert(power_tower(twos, 1000000000000LL) == power_tower(vector<long long>(twos.begin(), twos.begin() + 60), 1000000000000LL));

    /// A zero the evaluation reads must abort, one it never reaches is not validated
    for (auto bad : vector<vector<long long>>{{0}, {2, 0}, {2, 0, 3}, {3, 3, 3, 0, 5}}) assert(aborts(bad, 1000));

    PowerTower tower(1000);
    assert(tower.query({2, 3, 0}, 0, 1) == 8);
    assert(tower.query({7, 1, 0}, 0, 2) == 7);

    return 0;
}
