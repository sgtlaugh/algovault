#include "common.h"

#define main library_main
#include "../code_library/sparse_table.cpp"
#undef main

/// Range minimum queries against a direct scan, every (l, r) for small arrays and random ranges for large ones
int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it + 1 : stress::rand_int(1, 3000);
        vector<int> a(n);
        for (auto& x : a) x = stress::rand_int(INT_MIN, INT_MAX);
        SparseTable<int> table(a);

        if (n <= 64){
            for (int l = 0; l < n; l++){
                int expected = INT_MAX;
                for (int r = l; r < n; r++) expected = min(expected, a[r]), assert(table.query(l, r) == expected);
            }
        }
        for (int q = 0; q < 200; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            assert(table.query(l, r) == *min_element(a.begin() + l, a.begin() + r + 1));
        }
    }
    return 0;
}
