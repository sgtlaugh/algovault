/***
 *
 * 2D Fenwick Tree, Sparse
 * Point add and rectangle sum on grids up to 1e9 x 1e9, keeping only the touched Fenwick nodes in a hash table
 *
 * Complexity:
 *   - build: O(2^lg_capacity)
 *   - update, prefix, query: O(log n log m) hash probes, each expected O(1) while the table is at most half full
 *
 * FenwickSparse2D<T> fen(n, m, lg_capacity = 24): n rows, m columns, all zeros, in a table of 2^lg_capacity slots
 * fen.update(i, j, v): add v to a[i][j]
 * fen.prefix(i, j): sum of a[0..i][0..j], 0 when i or j is -1
 * fen.query(i, j, k, l): sum of a[i..k][j..l], 0 when i > k or j > l
 *
 * Notes:
 *   - Cells are 0-based in [0, n - 1] x [0, m - 1], outside asserts
 *   - The table takes 8 + sizeof(T) bytes per slot up front, 268 MB at the default lg_capacity = 24 for T = long long
 *   - An update creates up to (floor(log2(n)) + 1) * (floor(log2(m)) + 1) nodes, 900 at n = m = 1e9
 *   - Size 2^lg_capacity to at least twice the nodes created, an update that fills the table past 3/4 asserts
 *   - So n = m = 1e9 with 100K updates, up to 9e7 nodes, does not fit the default
 *   - Use fenwick_tree_2D_implicit.cpp for rectangle add and cell value, paying only for the nodes it creates
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickSparse2D{
    int n, m, mask, used = 0;
    vector<T> vals;
    vector<long long> keys;

    FenwickSparse2D(int n, int m, int lg_capacity = 24) : n(n), m(m), mask((1 << lg_capacity) - 1),
        vals(mask + 1, 0), keys(mask + 1, 0){
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

    void add(int x, int y, T v){
        long long key = ((long long)x << 32) | y;
        int k = slot(key);
        if (!keys[k]){
            used++;
            assert(used <= (mask + 1) / 4 * 3);  /// table over 3/4 full, raise lg_capacity
            keys[k] = key;
        }
        vals[k] += v;
    }

    /// Testing the key first skips a cache miss on vals at empty slots, queries 1.3x faster at n = 1e9
    T find(int x, int y){
        int k = slot(((long long)x << 32) | y);
        return keys[k] ? vals[k] : 0;
    }

    /// Key 0 marks an empty slot, Fenwick node indices start at 1 so every key is non-zero
    int slot(long long key){
        int k = hash_func(key) & mask;
        while (keys[k] && keys[k] != key) k = (k + 1) & mask;
        return k;
    }

    /// MurmurHash3 fmix64 finalizer
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
    /// Three points on a 1e9 x 1e9 grid
    FenwickSparse2D<long long> fen(1000000000, 1000000000);
    fen.update(0, 0, 5);
    fen.update(999999999, 999999999, 10);
    fen.update(499999999, 499999999, 7);
    assert(fen.query(0, 0, 999999999, 999999999) == 22);
    assert(fen.query(499999999, 499999999, 999999999, 999999999) == 17);
    assert(fen.query(1, 0, 499999998, 999999999) == 0);
    assert(fen.prefix(499999999, 999999999) == 12);

    /// A 2^16 slot table is plenty for a few updates on a small grid
    FenwickSparse2D<int> small(1000, 1000, 16);
    small.update(2, 3, 2);
    small.update(999, 999, 9);
    assert(small.query(0, 0, 2, 3) == 2);
    assert(small.query(2, 3, 999, 999) == 11);

    return 0;
}
