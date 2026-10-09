#include "common.h"

#define main library_main
#include "../code_library/stirling_numbers.cpp"
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

    /// Identities on a long row: sum of c(n, k) is n!, c(n, 1) = (n - 1)!, c(n, n - 1) = n (n - 1) / 2
    for (long long m : {1000000007LL, 1LL << 30}){
        int n = 100000;
        auto row = stirling_first(n, m);
        long long sum = 0, fact = 1, fact_minus_one = 1;
        for (long long x : row) sum = (sum + x) % m;
        for (int i = 1; i <= n; i++){
            fact = fact * i % m;
            if (i < n) fact_minus_one = fact_minus_one * i % m;
        }
        assert(sum == fact && row[1] == fact_minus_one && row[n - 1] == (long long)n * (n - 1) / 2 % m && row[n] == 1 % m);
    }
    return 0;
}
