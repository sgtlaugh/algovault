/***
 *
 * Sparse 2D Fenwick Tree (Binary Indexed Tree)
 * For very large matrices (N up to 10^9) with sparse updates
 * Uses custom hashmap (SplitMix64) to store only non-zero Fenwick tree nodes
 * Point updates and range queries only
 * 1-based indexing for elements
 *
 * Time: O(log N * log M) per update/query
 * Space: a fixed table of 2^lg_capacity slots, (8 + sizeof(T)) bytes each, allocated up front
 *   - The default lg_capacity = 24 is 16.7M slots, 268 MB for T = long long
 *   - Each update creates up to (log2(N) + 1) * (log2(M) + 1) nodes, 961 at N = M = 1e9
 *   - The table must hold every node created, a full table makes probing loop forever
 *   - Size it so 2^lg_capacity is at least twice Q * (log2(N) + 1) * (log2(M) + 1) to keep probes short
 *   - So N = M = 1e9 with 100K updates (up to 9e7 nodes) does not fit the default, it hangs
 *
 * Benchmarks (Codeforces server):
 *   - N = 100K, Q = 100K random operations: 0.8 seconds
 *   - N = 200K, Q = 200K random operations: 2.0 seconds
 *
 * vs implicit segment tree (fenwick_tree_2D_implicit.cpp):
 *   - Faster for very sparse updates (Q < 10K) with random access
 *   - The table costs its full size regardless of Q
 *   - Implicit better for moderate Q (10K-100K) due to cache locality
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickSparse2D{
    int n, m, mask;
    vector<T> tree;
    vector<long long> hashmap;

    FenwickSparse2D(int n = 0, int m = 0, int lg_capacity = 24) : n(n), m(m), mask((1 << lg_capacity) - 1),
        tree(mask + 1, 0), hashmap(mask + 1, 0){
    }

    inline unsigned long long hash_func(unsigned long long h){
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return h;
    }

    /// Key 0 marks an empty slot, i >= 1 keeps every key non-zero
    inline void add(int i, int j, T v){
        long long h = ((long long)i << 32) | j;
        int k = hash_func(h) & mask;
        while (hashmap[k] && hashmap[k] != h) k = (k + 1) & mask;
        hashmap[k] = h;
        tree[k] += v;
    }

    inline T find(int i, int j){
        long long h = ((long long)i << 32) | j;
        int k = hash_func(h) & mask;
        while (hashmap[k] && hashmap[k] != h) k = (k + 1) & mask;
        return (hashmap[k] ? tree[k] : 0);
    }

    void update(int i, int j, T v){
        for (int x = i; x <= n; x += x & -x){
            for (int y = j; y <= m; y += y & -y){
                add(x, y, v);
            }
        }
    }

    T query(int i, int j){
        if (i < 0 || j < 0 || i > n || j > m) return 0;

        T res = 0;
        for (int x = i; x > 0; x -= x & -x){
            for (int y = j; y > 0; y -= y & -y){
                res += find(x, y);
            }
        }
        return res;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        return query(k, l) - query(i - 1, l) - query(k, j - 1) + query(i - 1, j - 1);
    }
};

int main(){
    /// Sparse 2D Fenwick for large matrix
    auto fen = FenwickSparse2D<long long>(1000000000, 1000000000);

    fen.update(1, 1, 5);
    fen.update(1000000000, 1000000000, 10);
    fen.update(500000000, 500000000, 7);

    assert(fen.query(1, 1, 1, 1) == 5);
    assert(fen.query(1, 1, 1000000000, 1000000000) == 22);
    assert(fen.query(500000000, 500000000, 1000000000, 1000000000) == 17);

    /// A 2^16 slot table is plenty for a few updates on a small grid
    auto small = FenwickSparse2D<int>(1000, 1000, 16);
    small.update(3, 4, 2);
    small.update(1000, 1000, 9);

    assert(small.query(1, 1, 3, 4) == 2);
    assert(small.query(3, 4, 1000, 1000) == 11);

    return 0;
}
