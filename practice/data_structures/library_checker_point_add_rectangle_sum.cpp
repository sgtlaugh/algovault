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
 * Cells are 0-based, a cell outside [0, n - 1] x [0, m - 1] asserts
 * The table takes 8 + sizeof(T) bytes per slot up front, 268 MB at the default lg_capacity = 24 for T = long long
 * An update creates up to (floor(log2(n)) + 1) * (floor(log2(m)) + 1) nodes, 900 at n = m = 1e9
 * Size 2^lg_capacity to at least twice the nodes created to keep probes short
 * An update that fills the table past 3/4 asserts
 * So n = m = 1e9 with 100K updates, up to 9e7 nodes, does not fit the default
 * fenwick_tree_2D_implicit.cpp does rectangle add and cell value instead, paying only for the nodes it creates
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

    /// Most probes miss, testing the key first skips a cache miss on vals
    T find(int x, int y){
        int k = slot(((long long)x << 32) | y);
        return keys[k] ? vals[k] : 0;
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

    /// Key 0 marks an empty slot, Fenwick node indices start at 1 so every key is non-zero
    int slot(long long key){
        int k = hash_func(key) & mask;
        while (keys[k] && keys[k] != key) k = (k + 1) & mask;
        return k;
    }
};

int main(){
    int n, q;
    scanf("%d %d", &n, &q);

    vector<array<int, 3>> points(n);
    for (auto& p: points) scanf("%d %d %d", &p[0], &p[1], &p[2]);

    vector<array<int, 5>> queries(q);
    vector<int> xs, ys;
    for (auto& p: points) xs.push_back(p[0]), ys.push_back(p[1]);
    for (auto& qr: queries){
        scanf("%d", &qr[0]);
        if (qr[0] == 0){
            scanf("%d %d %d", &qr[1], &qr[2], &qr[3]);
            xs.push_back(qr[1]), ys.push_back(qr[2]);
        }
        else scanf("%d %d %d %d", &qr[1], &qr[2], &qr[3], &qr[4]);
    }

    /// Compressing to the at most 2e5 point coordinates cuts an update from 900 nodes to about 300,
    /// so the worst inputs create about 1.1e7 nodes, a third of the 2^25 slot table (512 MB)
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    auto rank_of = [](const vector<int>& v, int c){
        return int(lower_bound(v.begin(), v.end(), c) - v.begin());
    };

    FenwickSparse2D<long long> fen(xs.size(), ys.size(), 25);
    for (auto& p: points) fen.update(rank_of(xs, p[0]), rank_of(ys, p[1]), p[2]);
    for (auto& qr: queries){
        if (qr[0] == 0) fen.update(rank_of(xs, qr[1]), rank_of(ys, qr[2]), qr[3]);
        else{
            printf("%lld\n", fen.query(rank_of(xs, qr[1]), rank_of(ys, qr[2]), rank_of(xs, qr[3]) - 1, rank_of(ys, qr[4]) - 1));
        }
    }
    return 0;
}
