#include "../common.h"

#define main library_main
#include "../../code_library/misc/permutation_index.cpp"
#undef main

/// Lexicographic rank by counting smaller unused values at each position
long long brute_rank(const vector<int>& p){
    int n = p.size();
    long long res = 0;
    vector<long long> factorial(n + 1, 1);
    for (int i = 1; i <= n; i++) factorial[i] = factorial[i - 1] * i;
    for (int i = 0; i < n; i++){
        int smaller = 0;
        for (int j = i + 1; j < n; j++) smaller += p[j] < p[i];
        res += smaller * factorial[n - 1 - i];
    }
    return res;
}

int main(){
    /// Every permutation in order for small n: the index must count 0, 1, 2, ...
    for (int n = 1; n <= 8; n++){
        PermutationIndex idx(n);
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        long long expected = 0;
        do assert(idx.index(p) == expected++);
        while (next_permutation(p.begin(), p.end()));
    }

    for (int n = 9; n <= 12; n++){
        PermutationIndex idx(n);
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        for (long long it = 0; it < stress::scaled(50000); it++){
            shuffle(p.begin(), p.end(), stress::rng());
            assert(idx.index(p) == brute_rank(p));
        }
    }
    return 0;
}
