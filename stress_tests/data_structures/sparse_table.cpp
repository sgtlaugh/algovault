#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sparse_table.cpp"
#undef main

/// Range minimum queries against a direct scan, every (l, r) for small arrays and random ranges for large ones
int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it + 1 : stress::rand_int(1, 3000);
        vector<int> a(n);
        for (auto& x : a) x = stress::rand_int(INT_MIN, INT_MAX);
        if (it % 3 == 0) for (auto& x : a) x %= 4;
        SparseTable<int> table(a);
        LinearSparseTable<int> linear(a);

        if (n <= 200){
            for (int l = 0; l < n; l++){
                int expected = INT_MAX;
                for (int r = l; r < n; r++){
                    expected = min(expected, a[r]);
                    assert(table.query(l, r) == expected && linear.query(l, r) == expected);
                }
            }
        }
        for (int q = 0; q < 200; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            int expected = *min_element(a.begin() + l, a.begin() + r + 1);
            assert(table.query(l, r) == expected && linear.query(l, r) == expected);
        }
    }

    /// The linear table on its own at a size where every block boundary case appears many times
    vector<long long> big(1000000);
    for (auto& x : big) x = stress::rand_int(-1000000000000LL, 1000000000000LL);
    LinearSparseTable<long long> linear(big);
    for (int q = 0; q < 2000; q++){
        int l = stress::rand_int(0, 999999), r = stress::rand_int(l, min(999999, l + (q % 2 ? 200 : 1000000)));
        assert(linear.query(l, r) == *min_element(big.begin() + l, big.begin() + r + 1));
    }
    return 0;
}
