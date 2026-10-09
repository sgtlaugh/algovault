/***
 *
 * Chromatic Number
 * Minimum number of colors in a proper vertex coloring of an undirected graph
 *
 * Complexity: O(2^n n) time, 12 * 2^n bytes of memory (about 200 MB at n = 24)
 *
 * ChromaticNumber g(n): vertices 0..n - 1, n <= 24
 * g.add_edge(u, v): undirected edge, u != v, duplicates are fine
 * g.solve(): the chromatic number, 0 for the empty graph
 *
 * ind[S] = number of independent subsets of S. By inclusion-exclusion, the number of k-tuples of
 * independent sets whose union is V is sum over S of (-1)^(n - |S|) ind[S]^k, and it is nonzero iff
 * the graph is k-colorable. The sums are taken mod the prime 2^61 - 1: the answer is never too small,
 * and too large only if the true count at the chromatic number is a multiple of that prime
 *
 * Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ChromaticNumber{
    static constexpr uint64_t MOD = (1ULL << 61) - 1;

    int n;
    vector<int> adj;

    ChromaticNumber(int n) : n(n), adj(n, 0){
        assert(0 <= n && n <= 24);
    }

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        adj[u] |= 1 << v;
        adj[v] |= 1 << u;
    }

    int solve() const{
        if (n == 0) return 0;
        int full = 1 << n;

        vector<uint32_t> ind(full);
        ind[0] = 1;
        for (int mask = 1; mask < full; mask++){
            int u = __builtin_ctz(mask), rest = mask ^ (1 << u);
            ind[mask] = ind[rest] + ind[rest & ~adj[u]];
        }

        vector<uint64_t> power(full, 1);
        for (int k = 1; k < n; k++){
            uint64_t by_parity[2] = {0, 0};
            for (int mask = 0; mask < full; mask++){
                power[mask] = mul(power[mask], ind[mask]);
                uint64_t& sum = by_parity[__builtin_popcount(mask) & 1];
                sum += power[mask];
                if (sum >= MOD) sum -= MOD;
            }
            if (by_parity[0] != by_parity[1]) return k;
        }
        return n;
    }

private:
    static uint64_t mul(uint64_t a, uint64_t b){
        unsigned __int128 product = (unsigned __int128)a * b;
        uint64_t res = (uint64_t)(product & MOD) + (uint64_t)(product >> 61);
        return res >= MOD ? res - MOD : res;
    }
};

int main(){
    auto chromatic = [](int n, const vector<pair<int, int>>& edges){
        ChromaticNumber g(n);
        for (auto [u, v] : edges) g.add_edge(u, v);
        return g.solve();
    };

    assert(ChromaticNumber(5).solve() == 1);  /// no edges
    assert(chromatic(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}) == 2);
    assert(chromatic(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}) == 3);  /// an odd cycle
    assert(chromatic(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}) == 4);

    /// Grotzsch graph (Mycielskian of C5): no triangle, yet 4 colors, so cliques alone do not decide it
    vector<pair<int, int>> grotzsch;
    for (int i = 0; i < 5; i++){
        int prev = (i + 4) % 5, next = (i + 1) % 5;
        grotzsch.push_back({i, next});
        grotzsch.push_back({i + 5, prev});
        grotzsch.push_back({i + 5, next});
        grotzsch.push_back({i + 5, 10});
    }
    assert(chromatic(11, grotzsch) == 4);
    return 0;
}
