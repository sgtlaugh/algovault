/***
 *
 * 2D Fenwick Tree with Implicit Segment Trees for large matrices with sparse updates
 * Each Fenwick row uses implicit segment tree for columns - allocates nodes on demand
 * Range updates and point queries only
 * 0-based cells, rectangles [x1, x2] x [y1, y2] inclusive, a cell outside [0, n - 1] x [0, m - 1] asserts
 *
 * FenwickImplicit2D<long long> f(n, m): n rows, m columns, all zeros
 * f.update(x1, y1, x2, y2, v): add v to every cell of the rectangle, x1 > x2 or y1 > y2 is empty and does nothing
 * f.query(x, y): the value at cell (x, y)
 *
 * Complexity:
 *   - O(n) to construct, for the row roots
 *   - O(log n log m) per update/query, a range update is 4 point updates of O(log n) rows each creating up to
 *     log m + 1 nodes
 *   - O(q log n log m) nodes for q updates
 *
 * Nodes come from a pool of MAXNODES (16 bytes each for long long, ~205 MB by default), each range update takes up to
 * 4 * log(n) * log(m)
 * The pool is reserved up front but only memory for nodes actually created is touched
 * The row roots take another 4 * (n + 1) bytes, so the row count n is bounded by memory, m only by the pool
 *
 * vs sparse hashmap (fenwick_tree_2D_sparse.cpp):
 *   - Better for Q < 100K: no fixed hash table cost (268 MB by default)
 *   - Better cache locality: tree structure vs random hash positions
 *   - No hash collisions
 *   - Hashmap faster for very sparse (Q < 10K) with random access patterns
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickImplicit2D{
    static const int MAXNODES = 12800000;

    struct Node{
        T sum;
        int l, r;
    };

    int n, m;
    vector<int> root;
    vector<Node> nodes;

    /// reserve without filling, so pool pages are only touched once nodes are used, node 0 is the empty child
    FenwickImplicit2D(int n, int m) : n(n), m(m), root(n + 1, 0){
        nodes.reserve(MAXNODES);
        nodes.push_back({0, 0, 0});
    }

    void update(int x1, int y1, int x2, int y2, T v){
        if (x1 > x2 || y1 > y2) return;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m);
        update_fen(x1, y1, v);
        update_fen(x2 + 1, y1, -v);
        update_fen(x2 + 1, y2 + 1, v);
        update_fen(x1, y2 + 1, -v);
    }

    T query(int x, int y){
        assert(0 <= x && x < n && 0 <= y && y < m);
        T res = 0;
        for (int i = x + 1; i > 0; i -= i & -i){
            res += query_seg(root[i], 0, m - 1, y);
        }
        return res;
    }

private:
    int update_seg(int cur, int a, int b, int p, T v){
        if (!cur){
            assert((int)nodes.size() < MAXNODES);  /// out of pool nodes, raise MAXNODES
            cur = nodes.size();
            nodes.push_back({0, 0, 0});
        }
        nodes[cur].sum += v;

        if (a != b){
            int mid = a + (b - a) / 2;
            if (p <= mid) nodes[cur].l = update_seg(nodes[cur].l, a, mid, p, v);
            else nodes[cur].r = update_seg(nodes[cur].r, mid + 1, b, p, v);
        }
        return cur;
    }

    T query_seg(int cur, int a, int b, int r){
        if (!cur) return 0;
        if (r >= b) return nodes[cur].sum;

        int mid = a + (b - a) / 2;
        if (r <= mid) return query_seg(nodes[cur].l, a, mid, r);
        return nodes[nodes[cur].l].sum + query_seg(nodes[cur].r, mid + 1, b, r);
    }

    void update_fen(int x, int y, T v){
        if (y >= m) return;  /// column m from a range update ending at m - 1, otherwise it lands on leaf m - 1
        for (int i = x + 1; i <= n; i += i & -i){
            root[i] = update_seg(root[i], 0, m - 1, y, v);
        }
    }
};

int main(){
    /// Test with moderate size
    auto fen = FenwickImplicit2D<long long>(1000, 1000);

    fen.update(2, 4, 5, 6, 10);
    fen.update(3, 5, 7, 8, 5);

    assert(fen.query(2, 4) == 10);
    assert(fen.query(4, 5) == 15);
    assert(fen.query(6, 7) == 5);
    assert(fen.query(1, 3) == 0);

    /// Test range updates
    auto fen2 = FenwickImplicit2D<int>(100, 100);
    fen2.update(9, 19, 29, 39, 7);
    fen2.update(49, 59, 69, 79, 3);
    fen2.update(90, 95, 99, 99, 1);

    assert(fen2.query(19, 29) == 7);
    assert(fen2.query(59, 69) == 3);
    assert(fen2.query(0, 0) == 0);
    assert(fen2.query(34, 44) == 0);
    assert(fen2.query(99, 99) == 1);
    assert(fen2.query(89, 99) == 0);

    /// Test large coordinates
    auto fen3 = FenwickImplicit2D<long long>(100000, 100000);
    fen3.update(9999, 19999, 49999, 59999, 2);

    assert(fen3.query(29999, 39999) == 2);
    assert(fen3.query(9998, 19998) == 0);
    assert(fen3.query(49999, 59999) == 2);

    /// Rows and columns are sized separately, 2e9 columns cost no more row roots than 1e5 and a + b would overflow
    auto wide = FenwickImplicit2D<long long>(100000, 2000000000);
    wide.update(4, 6, 99999, 1999999999, 3);

    assert(wide.query(4, 6) == 3);
    assert(wide.query(99999, 1999999999) == 3);
    assert(wide.query(3, 6) == 0);
    assert(wide.query(4, 5) == 0);

    return 0;
}
