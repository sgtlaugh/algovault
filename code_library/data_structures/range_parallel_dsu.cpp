/***
 *
 * Range Parallel Union-Find
 * Union-find that joins two ranges position by position, a + i with b + i for every i < k, in one call
 *
 * Complexity: O((n log n + q) α(n)) for q connects in total, O(n log n) memory
 *
 * RangeParallelDSU dsu(n): nodes 0..n-1, each in its own set
 * dsu.connect(a, b, k): join a + i and b + i for 0 <= i < k, returns the number of merges, a + k <= n and b + k <= n
 * dsu.connect(a, b, k, on_merge): the same, calling on_merge(root, absorbed) once per merge of two sets
 * dsu.find(a): the root of a's set, the same node for every member
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 *
 * Notes:
 *   - Online: every connect takes effect at once, so find, same and size are exact between calls
 *   - on_merge runs right after each merge, with the new root and the absorbed set's old root
 *   - Per-set sums: on_merge(root, absorbed) adds sum[absorbed] into sum[root]
 *   - The ranges may overlap: connect(0, 1, n - 1) joins every node
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY dsu from code_library/data_structures/disjoint_set.cpp
struct DSU{
    vector<int> sz, parent;

    DSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        return true;
    }

    int find(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }

        return parent[i];
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }
};
/// END COPY dsu

struct RangeParallelDSU{
    int n;
    vector<DSU> layer;  /// layer[j] joins start a with start b: a + i ~ b + i for every i < 2^j

    RangeParallelDSU(int n) : n(n){
        for (int j = 0; (1 << j) <= n; j++) layer.emplace_back(n - (1 << j));
    }

    int connect(int a, int b, int k){
        return connect(a, b, k, [](int, int) {});
    }

    template <typename F>
    int connect(int a, int b, int k, F on_merge){
        assert(0 <= k && 0 <= a && 0 <= b && a + k <= n && b + k <= n);
        if (k == 0) return 0;

        /// Two blocks of length 2^j, the largest power of two <= k, cover [0, k) between them
        int j = __lg(k), merges = 0;
        link(j, a, b, on_merge, merges);
        link(j, a + k - (1 << j), b + k - (1 << j), on_merge, merges);
        return merges;
    }

    int find(int a){
        return layer[0].find(a);
    }

    bool same(int a, int b){
        return layer[0].same(a, b);
    }

    int size(int a){
        return layer[0].size(a);
    }

    /// A join that merges nothing is already implied, so only a real merge is pushed down: at most n per layer
    template <typename F>
    void link(int j, int a, int b, F& on_merge, int& merges){
        if (j == 0){
            int ra = layer[0].find(a), rb = layer[0].find(b);
            if (ra == rb) return;

            layer[0].connect(ra, rb);
            int root = layer[0].find(ra);
            on_merge(root, root == ra ? rb : ra);
            merges++;
            return;
        }

        if (!layer[j].connect(a, b)) return;
        int half = 1 << (j - 1);
        link(j - 1, a, b, on_merge, merges);
        link(j - 1, a + half, b + half, on_merge, merges);
    }
};

int main(){
    /// Six nodes: (0, 3), (1, 4) and (2, 5) pair up, then 0 ~ 1 and 1 ~ 2 join the three pairs
    RangeParallelDSU dsu(6);
    assert(dsu.connect(0, 3, 3) == 3);
    assert(dsu.same(1, 4));
    assert(!dsu.same(0, 1));
    assert(dsu.connect(0, 1, 2) == 2);
    assert(dsu.size(5) == 6);
    assert(dsu.connect(1, 3, 2) == 0);

    /// Weights 1..6 summed per set: 0 ~ 2 and 1 ~ 3 give sets {0, 2} of weight 4 and {1, 3} of weight 6
    vector<long long> weight = {1, 2, 3, 4, 5, 6};
    RangeParallelDSU sums(6);
    sums.connect(0, 2, 2, [&](int root, int absorbed){ weight[root] += weight[absorbed]; });
    assert(weight[sums.find(0)] == 4);
    assert(weight[sums.find(3)] == 6);
    assert(weight[sums.find(4)] == 5);

    return 0;
}
