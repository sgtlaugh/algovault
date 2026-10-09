#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/disjoint_sparse_table.cpp"
#undef main

/// combine is x + y, so strings make it concatenation: associative but not commutative, operand order must be exact
int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 40 : 600);
        if (it % 7 == 0) n = 1 << stress::rand_int(0, 9);  /// block boundaries line up with n

        vector<string> words(n);
        vector<long long> nums(n);
        for (int i = 0; i < n; i++){
            words[i] = string(1, 'a' + i % 26) + (i % 3 ? "" : to_string(i));
            nums[i] = stress::rand_int(-1000000000000LL, 1000000000000LL);
        }

        DisjointST<string> concat(words, "");
        DisjointST<long long> sum(nums, 0);
        for (int q = 0; q < 50; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, min(n - 1, l + (q % 2 ? 5 : n)));
            string expected;
            for (int i = l; i <= r; i++) expected += words[i];
            assert(concat.query(l, r) == expected);
            assert(sum.query(l, r) == accumulate(nums.begin() + l, nums.begin() + r + 1, 0LL));
        }
    }

    return 0;
}
