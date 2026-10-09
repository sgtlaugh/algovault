#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/partition_numbers.cpp"
#undef main

/// p(0..n) modulo m by the O(n^2) coin DP over part sizes 1..n
vector<long long> coin_dp(int n, long long m){
    vector<unsigned long long> p(n + 1, 0);
    p[0] = 1 % m;
    for (int c = 1; c <= n; c++){
        for (int i = c; i <= n; i++){
            p[i] += p[i - c];
            if (p[i] >= (unsigned long long)m) p[i] -= m;
        }
    }
    return vector<long long>(p.begin(), p.end());
}

/// Compares partition_numbers against the coin DP, and Ramanujan's congruences at large n
int main(){
    for (long long m = 1; m <= 40; m++){
        auto expected = coin_dp(80, m);
        for (int n = 0; n <= 80; n++){
            assert(partition_numbers(n, m) == vector<long long>(expected.begin(), expected.begin() + n + 1));
        }
    }

    vector<long long> moduli = {1, 2, 1000000007, 998244353, (1LL << 62) - 1, 1LL << 62, 999999999999999989LL};
    for (long long it = 0; it < stress::scaled(10); it++){
        moduli.push_back(stress::rand_int(1, 1LL << 62));
        moduli.push_back(stress::rand_int(1, 1000));
    }

    for (long long m : moduli){
        int n = stress::rand_int(1500, 3000);
        assert(partition_numbers(n, m) == coin_dp(n, m));
    }

    /// p(5k + 4) = 0 mod 5, p(7k + 5) = 0 mod 7, p(11k + 6) = 0 mod 11, all three read off one run modulo 385
    int n = 100000;
    auto p = partition_numbers(n, 385);
    for (int i = 4; i <= n; i += 5) assert(p[i] % 5 == 0);
    for (int i = 5; i <= n; i += 7) assert(p[i] % 7 == 0);
    for (int i = 6; i <= n; i += 11) assert(p[i] % 11 == 0);

    return 0;
}
