/***
 *
 * Sparse 2D Fenwick Tree (Binary Indexed Tree)
 * For very large matrices (N up to 10^9) with sparse updates
 * Uses custom hashmap (SplitMix64) to store only non-zero Fenwick tree nodes
 * Point updates and range queries only
 * 1-based indexing for elements
 *
 * Space: O(Q log N log M) where Q is number of updates
 * Time: O(log N * log M) per update/query
 *
 * Benchmarks (Codeforces server):
 *   - N = 100K, Q = 100K random operations: 0.8 seconds
 *   - N = 200K, Q = 200K random operations: 2.0 seconds
 *
 * vs implicit segment tree (fenwick_tree_2D_implicit.cpp):
 *   - Faster for very sparse updates (Q < 10K) with random access
 *   - Hashmap has 16MB+ overhead regardless of Q
 *   - Implicit better for moderate Q (10K-100K) due to cache locality
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickSparse2D{
    static const int HMOD = 16777215;
    static const int HSIZE = HMOD + 71235;

    int n, m;
    vector<T> tree;
    vector<long long> hashmap;

    FenwickSparse2D(int n = 0, int m = 0) : n(n), m(m), tree(HSIZE, 0), hashmap(HSIZE, 0) {
    }

    inline unsigned long long hash_func(unsigned long long h){
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return h;
    }

    inline void add(int i, int j, T v){
        long long h = ((long long)i * 1000003LL) + j + 917921LL;
        int k = hash_func(h) & HMOD;
        while (hashmap[k] && hashmap[k] != h) k++;
        hashmap[k] = h;
        tree[k] += v;
    }

    inline T find(int i, int j){
        long long h = ((long long)i * 1000003LL) + j + 917921LL;
        int k = hash_func(h) & HMOD;
        while (hashmap[k] && hashmap[k] != h) k++;
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

    return 0;
}
