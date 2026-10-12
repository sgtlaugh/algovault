/***
 *
 * Sparse 2D Fenwick Tree (Binary Indexed Tree)
 * For very large matrices (N up to 10^9) with sparse updates
 * Uses custom hashmap (MurmurHash3 fmix64 finalizer) to store only non-zero Fenwick tree nodes
 * Point updates and range queries only
 * 0-based cells, rectangles [i, k] x [j, l] inclusive, a cell outside [0, n - 1] x [0, m - 1] asserts
 * i > k or j > l is an empty rectangle and query returns 0
 *
 * FenwickSparse2D<long long> f(n, m, lg_capacity);
 * f.update(i, j, v): add v to cell [i, j]
 * f.prefix(i, j): sum of rectangle [0, 0] to [i, j], 0 when i or j is -1
 * f.query(i, j, k, l): sum of rectangle [i, j] to [k, l]
 *
 * Complexity: O(2^lg_capacity) to construct, O(log N log M) hash probes per update/query, each expected O(1) while the table is at most half full
 * Space: a fixed table of 2^lg_capacity slots, (8 + sizeof(T)) bytes each, allocated up front
 *   - The default lg_capacity = 24 is 16.7M slots, 268 MB for T = long long
 *   - Each update creates up to (floor(log2(N)) + 1) * (floor(log2(M)) + 1) nodes, 900 at N = M = 1e9
 *   - The table must hold every node created, an update that fills it past 3/4 asserts, raise lg_capacity
 *   - Size it so 2^lg_capacity is at least twice Q * (log2(N) + 1) * (log2(M) + 1) to keep probes short
 *   - So N = M = 1e9 with 100K updates (up to 9e7 nodes) does not fit the default, it asserts
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
    int n, m, mask, used = 0;
    vector<T> tree;
    vector<long long> hashmap;

    FenwickSparse2D(int n, int m, int lg_capacity = 24) : n(n), m(m), mask((1 << lg_capacity) - 1),
        tree(mask + 1, 0), hashmap(mask + 1, 0){
    }

    void update(int i, int j, T v){
        assert(0 <= i && i < n && 0 <= j && j < m);
        for (int x = i + 1; x <= n; x += x & -x){
            for (int y = j + 1; y <= m; y += y & -y){
                add(x, y, v);
            }
        }
    }

    T prefix(int i, int j){
        assert(-1 <= i && i < n && -1 <= j && j < m);
        T res = 0;
        for (int x = i + 1; x > 0; x -= x & -x){
            for (int y = j + 1; y > 0; y -= y & -y){
                res += find(x, y);
            }
        }
        return res;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        assert(0 <= i && k < n && 0 <= j && l < m);
        return prefix(k, l) - prefix(i - 1, l) - prefix(k, j - 1) + prefix(i - 1, j - 1);
    }

private:
    /// Key 0 marks an empty slot, Fenwick node indices start at 1 so every key is non-zero
    void add(int x, int y, T v){
        long long h = ((long long)x << 32) | y;
        int k = hash_func(h) & mask;
        while (hashmap[k] && hashmap[k] != h) k = (k + 1) & mask;
        if (!hashmap[k]){
            used++;
            assert(used <= (mask + 1) / 4 * 3);  /// table over 3/4 full, raise lg_capacity
            hashmap[k] = h;
        }
        tree[k] += v;
    }

    T find(int x, int y){
        long long h = ((long long)x << 32) | y;
        int k = hash_func(h) & mask;
        while (hashmap[k] && hashmap[k] != h) k = (k + 1) & mask;
        return (hashmap[k] ? tree[k] : 0);
    }

    unsigned long long hash_func(unsigned long long h){
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return h;
    }
};

int main(){
    /// Sparse 2D Fenwick for large matrix
    auto fen = FenwickSparse2D<long long>(1000000000, 1000000000);

    fen.update(0, 0, 5);
    fen.update(999999999, 999999999, 10);
    fen.update(499999999, 499999999, 7);

    assert(fen.query(0, 0, 0, 0) == 5);
    assert(fen.query(0, 0, 999999999, 999999999) == 22);
    assert(fen.query(499999999, 499999999, 999999999, 999999999) == 17);
    assert(fen.query(1, 0, 499999998, 999999999) == 0);
    assert(fen.prefix(499999999, 999999999) == 12);
    assert(fen.prefix(-1, 999999999) == 0);

    /// A 2^16 slot table is plenty for a few updates on a small grid
    auto small = FenwickSparse2D<int>(1000, 1000, 16);
    small.update(2, 3, 2);
    small.update(999, 999, 9);

    assert(small.query(0, 0, 2, 3) == 2);
    assert(small.query(2, 3, 999, 999) == 11);
    assert(small.query(3, 3, 2, 3) == 0);

    return 0;
}
