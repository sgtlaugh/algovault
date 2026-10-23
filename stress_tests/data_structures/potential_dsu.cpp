#include "../common.h"

/// The library built with NDEBUG, as judges may do, so work done only inside its asserts would vanish
#define NDEBUG
#include <assert.h>
/// Under NDEBUG the values diff and main() compute only for asserts go unused, a warning judges do not treat as an error
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#define main library_main
#include "../../code_library/data_structures/potential_dsu.cpp"
#undef main
#pragma GCC diagnostic pop
#undef NDEBUG
#include <assert.h>

/// Union by size keeps every node within log2(component size) parent hops of its root
template <typename Dsu>
void check_depths(const Dsu& dsu){
    for (int i = 0; i < (int)dsu.parent.size(); i++){
        int depth = 0, root = i;
        while (dsu.parent[root] != root) root = dsu.parent[root], depth++;
        assert(depth < 31 && (1 << depth) <= dsu.sz[root]);
    }
}

void relabel(vector<int>& label, int a, int b){
    int from = label[a], to = label[b];
    if (from != to) for (auto& x : label) if (x == from) x = to;
}

/// Relations generated from a hidden assignment, some perturbed into contradictions once the pair is connected
int main(){
    for (int n : {1000, 100000}){
        PotentialDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.relate(i, i + 1, 1));
        check_depths(chain);
        assert(chain.diff(0, n) == n && chain.diff(n, 0) == -n && !chain.relate(n, 0, n));
    }

    /// diff must compress paths itself: 3 hangs two hops below the root and no other call finds it
    {
        PotentialDSU dsu(3);
        assert(dsu.relate(0, 1, 1));
        assert(dsu.relate(2, 3, 1));
        assert(dsu.relate(0, 2, 10));
        assert(dsu.diff(0, 3) == 11);
    }

    /// Pairwise merges of equal blocks never compress a path, so node i sits popcount(i) hops deep
    {
        const int n = (1 << 17) - 1;
        const long long mod = 1e18;
        vector<long long> x(n + 1);
        for (auto& v : x) v = stress::rand_int(0, mod - 1);

        PotentialDSU binomial(n, mod);
        for (int step = 1; step <= n; step *= 2){
            for (int i = 0; i + step <= n; i += 2 * step) assert(binomial.relate(i, i + step, x[i + step] - x[i]));
        }
        check_depths(binomial);

        /// The deepest node sums 17 potentials near 1e18, which overflows unless reduced per hop
        for (int i = n; i >= 0; i--) assert(binomial.diff(i, 0) == ((x[0] - x[i]) % mod + mod) % mod);
    }

    const long long lo = LLONG_MIN / 2, hi = LLONG_MAX / 2;
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 300);
        long long mod = vector<long long>{0, 0, 1, 2, 3, stress::rand_int(4, 1e9), (long long)1e18}[it % 7];
        bool extreme = mod == 0 && it % 2;

        vector<long long> x(n + 1);
        for (auto& v : x) v = mod ? stress::rand_int(0, mod - 1) : extreme ? stress::rand_int(lo, hi) : stress::rand_int(-1e18, 1e18);
        auto truth = [&](int a, int b){
            if (!mod) return x[b] - x[a];
            return ((x[b] - x[a]) % mod + mod) % mod;
        };

        PotentialDSU dsu(n, mod);
        vector<int> label(n + 1);
        iota(label.begin(), label.end(), 0);

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            long long d = x[b] - x[a];
            bool lie = label[a] == label[b] && mod != 1 && stress::rand_int(0, 2) == 0;

            if (lie && mod) d += stress::rand_int(1, mod - 1);
            if (lie && !mod) d += d >= 0 ? -stress::rand_int(1, 5) : stress::rand_int(1, 5);
            if (mod) d += mod * stress::rand_int(-3, 3);

            assert(dsu.relate(a, b, d) == !lie);
            relabel(label, a, b);

            int i = stress::rand_int(0, n), j = stress::rand_int(0, n);
            assert(dsu.same(i, j) == (label[i] == label[j]));
            assert(dsu.size(i) == count(label.begin(), label.end(), label[i]));
            if (label[i] == label[j]) assert(dsu.diff(i, j) == truth(i, j));
        }

        for (int i = 0; i <= n; i++){
            for (int j = 0; j <= n; j++) if (label[i] == label[j]) assert(dsu.diff(i, j) == truth(i, j));
        }
        check_depths(dsu);
    }

    return 0;
}
