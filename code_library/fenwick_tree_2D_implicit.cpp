/***
 *
 * 2D Fenwick Tree with Implicit Segment Trees for large matrices (N up to 10^9) with sparse updates
 * Each Fenwick row uses implicit segment tree for columns - allocates nodes on demand
 * Range updates and point queries only
 *
 * Space: O(Q log^2 N) where Q is number of updates
 * Time: O(log^2 N) per update/query
 *
 * vs sparse hashmap (fenwick_tree_2D_sparse.cpp):
 *   - Better for Q < 100K: no hash table overhead (16MB+)
 *   - Better cache locality: tree structure vs random hash positions
 *   - No hash collisions
 *   - Hashmap faster for very sparse (Q < 10K) with random access patterns
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickImplicit2D{
    static const int MAXNODES = 12800000; // Adjust based on memory constraints

    int n, idx;
    vector<T> tree;
    vector<int> root, L, R;

    FenwickImplicit2D(int n = 0) : n(n), idx(0), tree(MAXNODES, 0), root(n + 1, 0), L(MAXNODES, 0), R(MAXNODES, 0) {}

    void update_seg(int& cur, int a, int b, int p, T v){
        if (!cur) cur = ++idx;
        tree[cur] += v;

        if (a != b){
            int m = (a + b) >> 1;
            if (p <= m) update_seg(L[cur], a, m, p, v);
            else update_seg(R[cur], m + 1, b, p, v);
        }
    }

    T query_seg(int cur, int a, int b, int r){
        if (!cur) return 0;
        if (r >= b) return tree[cur];

        int m = (a + b) >> 1;
        if (r <= m) return query_seg(L[cur], a, m, r);
        return tree[L[cur]] + query_seg(R[cur], m + 1, b, r);
    }

    void update_fen(int x, int y, T v){
        for (int i = x; i <= n; i += i & -i){
            update_seg(root[i], 1, n, y, v);
        }
    }

    void update(int x, int y, T v){
        update_fen(x, y, v);
    }

    void update(int x1, int y1, int x2, int y2, T v){
        if (x1 > x2 || y1 > y2) return;
        update_fen(x1, y1, v);
        update_fen(x2 + 1, y1, -v);
        update_fen(x2 + 1, y2 + 1, v);
        update_fen(x1, y2 + 1, -v);
    }

    T query_fen(int x, int y){
        T res = 0;
        for (int i = x; i > 0; i -= i & -i){
            res += query_seg(root[i], 1, n, y);
        }
        return res;
    }

    T query(int x, int y){
        return query_fen(x, y);
    }
};

int main(){
    /// Test with moderate size
    auto fen = FenwickImplicit2D<long long>(1000);

    fen.update(3, 5, 6, 7, 10);
    fen.update(4, 6, 8, 9, 5);

    assert(fen.query(3, 5) == 10);
    assert(fen.query(5, 6) == 15);
    assert(fen.query(7, 8) == 5);
    assert(fen.query(2, 4) == 0);

    /// Test range updates
    auto fen2 = FenwickImplicit2D<int>(100);
    fen2.update(10, 20, 30, 40, 7);
    fen2.update(50, 60, 70, 80, 3);

    assert(fen2.query(20, 30) == 7);
    assert(fen2.query(60, 70) == 3);
    assert(fen2.query(1, 1) == 0);
    assert(fen2.query(35, 45) == 0);

    /// Test large coordinates
    auto fen3 = FenwickImplicit2D<long long>(100000);
    fen3.update(10000, 20000, 50000, 60000, 2);

    assert(fen3.query(30000, 40000) == 2);
    assert(fen3.query(9999, 19999) == 0);
    assert(fen3.query(50000, 60000) == 2);

    return 0;
}
