#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sparse_table.cpp"
#undef main

void check_rejects(int n, int l, int r){
    vector<int> a(n, 7);
    stress::expect_assert("data_structures/sparse_table.cpp", [&]{ SparseTable<int>(a).query(l, r); });
    stress::expect_assert("data_structures/sparse_table.cpp", [&]{ LinearSparseTable<int>(a).query(l, r); });
}

/// Range min, max and gcd queries against a direct scan, every (l, r) for small arrays and random ranges for large ones
int main(){
    check_rejects(9, 3, 1);
    check_rejects(9, -1, 2);
    check_rejects(9, 0, 9);
    check_rejects(200, 150, 200);  /// past the 64-wide blocks of LinearSparseTable
    check_rejects(0, 0, 0);

    auto max_op = [](int x, int y){ return max(x, y); };
    auto gcd_op = [](int x, int y){ return __gcd(x, y); };

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it + 1 : stress::rand_int(1, 3000);
        vector<int> a(n), pos(n);
        for (auto& x : a) x = stress::rand_int(INT_MIN, INT_MAX);
        if (it % 3 == 0) for (auto& x : a) x %= 4;
        for (int i = 0; i < n; i++) pos[i] = it % 2 ? stress::rand_int(0, 60) : abs(a[i] % 1000) * stress::rand_int(1, 12);

        SparseTable<int> table(a);
        LinearSparseTable<int> linear(a);
        SparseTable mx(a, max_op);
        LinearSparseTable<int, greater<int>> linear_mx(a);
        SparseTable g(pos, gcd_op);

        auto check = [&](int l, int r, int lo, int hi, int gc){
            assert(table.query(l, r) == lo && linear.query(l, r) == lo);
            assert(mx.query(l, r) == hi && linear_mx.query(l, r) == hi);
            assert(g.query(l, r) == gc);
        };

        if (n <= 200){
            for (int l = 0; l < n; l++){
                int lo = INT_MAX, hi = INT_MIN, gc = 0;
                for (int r = l; r < n; r++){
                    lo = min(lo, a[r]), hi = max(hi, a[r]), gc = gcd(gc, pos[r]);
                    check(l, r, lo, hi, gc);
                }
            }
        }

        for (int q = 0; q < 200; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            int gc = 0;
            for (int i = l; i <= r; i++) gc = gcd(gc, pos[i]);
            check(l, r, *min_element(a.begin() + l, a.begin() + r + 1), *max_element(a.begin() + l, a.begin() + r + 1), gc);
        }
    }

    /// The linear table on its own at a size where every block boundary case appears many times
    vector<long long> big(1000000);
    for (auto& x : big) x = stress::rand_int(-1000000000000LL, 1000000000000LL);
    LinearSparseTable<long long> linear(big);
    LinearSparseTable<long long, greater<long long>> linear_mx(big);

    for (int q = 0; q < 2000; q++){
        int l = stress::rand_int(0, 999999), r = stress::rand_int(l, min(999999, l + (q % 2 ? 200 : 1000000)));
        assert(linear.query(l, r) == *min_element(big.begin() + l, big.begin() + r + 1));
        assert(linear_mx.query(l, r) == *max_element(big.begin() + l, big.begin() + r + 1));
    }

    return 0;
}
