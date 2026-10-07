#include "common.h"

#define main library_main
#include "../code_library/permutation_rank.cpp"
#undef main

int main(){
    /// Every permutation in lexicographic order, ranks must be 1 .. n!
    for (int n = 0; n <= 7; n++){
        vector<int> p(n);
        iota(p.begin(), p.end(), 1);
        long long k = 1;
        do {
            assert(find_rank(p) == k);
            assert(find_permutation(n, k) == p);
            k++;
        } while (next_permutation(p.begin(), p.end()));
    }

    /// Up to 20! the ranks fit in 64 bits: round trip random ranks, the extremes included
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 20);
        long long total = 1;
        for (int i = 2; i <= n; i++) total *= i;
        long long k = it % 5 == 0 ? (it % 2 ? 1 : total) : stress::rand_int(1, total);

        auto p = find_permutation(n, k);
        auto sorted = p;
        sort(sorted.begin(), sorted.end());
        for (int i = 0; i < n; i++) assert(sorted[i] == i + 1);
        assert(find_rank(p) == k);
    }
    return 0;
}
