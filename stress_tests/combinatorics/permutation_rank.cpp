#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/permutation_rank.cpp"
#undef main

/// Factorial number system with vector erase, an independent O(n^2) unrank of 0..n-1
vector<int> naive_unrank(int n, __int128 k){
    vector<__int128> fact(n + 1, 1);
    for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i;
    vector<int> pool(n), res;
    iota(pool.begin(), pool.end(), 0);
    for (int i = n - 1; i >= 0; i--){
        int q = k / fact[i];
        k %= fact[i];
        res.push_back(pool[q]);
        pool.erase(pool.begin() + q);
    }
    return res;
}

/// Counts unused smaller values left to right in 128 bits, so 21 elements do not overflow
__int128 naive_rank(const vector<int>& p){
    int n = p.size();
    vector<bool> used(n, false);
    __int128 res = 0;
    for (int i = 0; i < n; i++){
        int smaller = 0;
        for (int v = 0; v < p[i]; v++) smaller += !used[v];
        used[p[i]] = true;
        res = res * (n - i) + smaller;
    }
    return res;
}

int main(){
    /// Every permutation in lexicographic order, ranks must be 0 .. n! - 1
    for (int n = 0; n <= 7; n++){
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        long long k = 0;
        do {
            assert(permutation_rank(p) == k);
            assert(kth_permutation(n, k) == p);
            k++;
        } while (next_permutation(p.begin(), p.end()));
    }

    /// Up to 20! the ranks fit in 64 bits: random ranks, the extremes included, on shifted distinct values
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 20);
        long long total = 1;
        for (int i = 2; i <= n; i++) total *= i;
        long long k = it % 5 == 0 ? (it % 2 ? 0 : total - 1) : stress::rand_int(0, total - 1);

        auto p = kth_permutation(n, k);
        assert(p == naive_unrank(n, k));
        assert(permutation_rank(p) == k);

        long long shift = stress::rand_int(-1000000000, 1000000000);
        vector<long long> shifted(p.begin(), p.end());
        for (auto& x: shifted) x = x * 3 + shift;
        assert(permutation_rank(shifted) == k);
    }

    /// n > 20: the prefix is the identity and the 21-element tail ranks back to k
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 100 ? 21 + it : stress::rand_int(21, it % 50 == 0 ? 100000 : 300);
        long long k = it % 7 == 0 ? LLONG_MAX - it % 3 : stress::rand_int(0, LLONG_MAX);

        auto p = kth_permutation(n, k);
        assert((int)p.size() == n);
        for (int i = 0; i < n - 21; i++) assert(p[i] == i);
        vector<int> tail(p.end() - 21, p.end());
        for (auto& x: tail) x -= n - 21;
        assert(tail == naive_unrank(21, k));
        assert(naive_rank(tail) == k);
    }

    return 0;
}
