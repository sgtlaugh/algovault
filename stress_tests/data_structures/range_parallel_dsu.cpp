#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/range_parallel_dsu.cpp"
#undef main

/// The reference: a plain DSU joined k times per call, with a set of current roots and per-set weight sums
struct Naive{
    DSU dsu;
    set<int> roots;
    vector<long long> weight;

    Naive(int n, const vector<long long>& w) : dsu(n), weight(w){
        for (int i = 0; i < n; i++) roots.insert(i);
    }

    int connect(int a, int b, int k){
        int merges = 0;
        for (int i = 0; i < k; i++){
            int ra = dsu.find(a + i), rb = dsu.find(b + i);
            if (ra == rb) continue;

            dsu.connect(ra, rb);
            int root = dsu.find(ra), gone = root == ra ? rb : ra;
            weight[root] += weight[gone], roots.erase(gone), merges++;
        }

        return merges;
    }
};

/// Every merge reported by on_merge is checked as it happens: both arguments were roots, and root now holds absorbed
void run(int n, int q, int max_k, bool all_pairs){
    vector<long long> w(n);
    for (auto& x : w) x = stress::rand_int(-1000000, 1000000);

    RangeParallelDSU dsu(n);
    Naive naive(n, w);
    set<int> roots;
    for (int i = 0; i < n; i++) roots.insert(i);

    for (int it = 0; it < q; it++){
        int k = stress::rand_int(0, min(n, max_k));
        int a = stress::rand_int(0, n - k), b = stress::rand_int(0, n - k);
        if (it % 5 == 0) b = min(n - k, a + (int)stress::rand_int(0, 3));  /// overlapping ranges, period 1 to 3

        int calls = 0;
        int merges = dsu.connect(a, b, k, [&](int root, int absorbed){
            assert(root != absorbed && roots.count(root) && roots.count(absorbed));
            assert(dsu.find(absorbed) == root && dsu.find(root) == root);
            roots.erase(absorbed), w[root] += w[absorbed], calls++;
        });
        assert(merges == calls);
        assert(merges == naive.connect(a, b, k));
        assert(roots == naive.roots);

        if (all_pairs || it % 50 == 0){
            for (int i = 0; i < n; i++){
                assert(dsu.size(i) == naive.dsu.size(i));
                assert(w[dsu.find(i)] == naive.weight[naive.dsu.find(i)]);
                for (int j = all_pairs ? 0 : i; j < (all_pairs ? n : min(n, i + 3)); j++){
                    assert(dsu.same(i, j) == naive.dsu.same(i, j));
                }
            }
        }
    }
}

int main(){
    stress::expect_assert("data_structures/range_parallel_dsu.cpp", []{ RangeParallelDSU dsu(5); dsu.connect(0, 3, 3); });
    stress::expect_assert("data_structures/range_parallel_dsu.cpp", []{ RangeParallelDSU dsu(5); dsu.connect(-1, 2, 2); });

    RangeParallelDSU empty(0);
    assert(empty.connect(0, 0, 0) == 0);

    /// Small n with every pair compared after each call, so a wrong push down to a lower layer shows at once
    for (int t = 0; t < stress::scaled(3000); t++){
        int n = stress::rand_int(1, 20);
        run(n, stress::rand_int(1, 30), n, true);
    }

    /// Medium n with long ranges, so joins start at the top layers and cross many block boundaries
    for (int t = 0; t < stress::scaled(30); t++){
        int n = stress::rand_int(100, 2000);
        run(n, 400, t % 2 ? n : 8, false);
    }

    /// Two hundred thousand nodes with long joins, the Library Checker scale, against the merge count bound
    int n = 200000, total = 0;
    RangeParallelDSU big(n);
    for (int it = 0; it < stress::scaled(200000); it++){
        int k = stress::rand_int(0, n);
        total += big.connect(stress::rand_int(0, n - k), stress::rand_int(0, n - k), k);
    }
    assert(total <= n - 1);
    int components = 0;
    for (int i = 0; i < n; i++) components += big.find(i) == i;
    assert(total == n - components);

    return 0;
}
