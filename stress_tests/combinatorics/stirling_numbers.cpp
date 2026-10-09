#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/combinatorics/stirling_numbers.cpp"
#undef main

bool is_prime(long long x){
    if (x < 2) return false;
    for (long long d = 2; d * d <= x; d++){
        if (x % d == 0) return false;
    }
    return true;
}

/// Rows 0..n of both kinds by their recurrences
void recurrence(int n, long long m, vector<vector<long long>>& first, vector<vector<long long>>& second){
    first.assign(n + 1, vector<long long>(n + 1, 0));
    second = first;
    first[0][0] = second[0][0] = 1 % m;
    for (int i = 0; i < n; i++){
        for (int k = 0; k <= i + 1; k++){
            first[i + 1][k] = ((k <= i ? (long long)i * first[i][k] : 0) + (k ? first[i][k - 1] : 0)) % m;
            second[i + 1][k] = ((k <= i ? (long long)k * second[i][k] : 0) + (k ? second[i][k - 1] : 0)) % m;
        }
    }
}

int main(){
    const long long BIG_PRIME = 1073741789;  /// largest prime below 2^30
    vector<long long> moduli = {1, 2, 6, 7, 1000000007, BIG_PRIME, 1LL << 30, 999999937};
    for (long long it = 0; it < stress::scaled(12); it++) moduli.push_back(stress::rand_int(1, 1LL << 30));

    for (long long m : moduli){
        int n = 120;
        vector<int> rows;
        for (int i = 0; i <= 30; i++) rows.push_back(i);
        for (int i : {60, 90, 120}) rows.push_back(i);

        vector<vector<long long>> first, second;
        recurrence(n, m, first, second);
        for (int i : rows){
            vector<long long> row_first(first[i].begin(), first[i].begin() + i + 1);
            assert(stirling_first(i, m) == row_first);
            if (is_prime(m) && i < m){
                vector<long long> row_second(second[i].begin(), second[i].begin() + i + 1);
                assert(stirling_second(i, m) == row_second);
            }
        }
    }

    for (long long m : {1000000007LL, BIG_PRIME}){
        int n = 1500;
        vector<vector<long long>> first, second;
        recurrence(n, m, first, second);
        assert(stirling_first(n, m) == first[n]);
        assert(stirling_second(n, m) == second[n]);
    }

    /// Long rows against their generating identities at a random point r, which catches any wrong coefficient:
    /// sum of c(n, k) r^k = r (r + 1) ... (r + n - 1) and sum of S(n, k) r (r - 1) ... (r - k + 1) = r^n
    {
        int n = 30000;
        long long m = 1000000007, r = stress::rand_int(1, m - 1), lhs = 0, rhs = 1, power = 1;
        auto row = stirling_first(n, m);
        for (long long x : row) lhs = (lhs + x * power) % m, power = power * r % m;
        for (int i = 0; i < n; i++) rhs = rhs * ((r + i) % m) % m;
        assert(lhs == rhs);
    }

    {
        int n = 120000;  /// the product length reaches 2^18
        long long p = BIG_PRIME, r = stress::rand_int(n + 1, p - 1), lhs = 0, rhs = 1, falling = 1;
        auto row = stirling_second(n, p);
        for (int k = 0; k <= n; k++) lhs = (lhs + row[k] * falling) % p, falling = falling * (r - k) % p;
        for (int i = 0; i < n; i++) rhs = rhs * r % p;
        assert(lhs == rhs);
    }

    /// A modulus below 2 must abort, stirling_second(0, 1) used to spin forever in pow_mod(1, -1)
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        alarm(5);
        assert(freopen("/dev/null", "w", stderr));
        stirling_second(0, 1);
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);

    return 0;
}
