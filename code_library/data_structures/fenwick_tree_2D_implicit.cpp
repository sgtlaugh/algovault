/***
 *
 * 2D Fenwick Tree, Implicit
 * Rectangle add and cell value on grids too large to allocate, each Fenwick row a segment tree over columns built on demand
 *
 * Complexity:
 *   - build: O(n) for the row roots
 *   - update, query: O(log n * log m)
 *
 * FenwickImplicit2D<T> fen(n, m): n rows, m columns, all zeros
 * fen.update(x1, y1, x2, y2, v): add v to a[x1..x2][y1..y2], nothing when x1 > x2 or y1 > y2
 * fen.query(x, y): a[x][y]
 *
 * Cells are 0-based, a cell outside [0, n - 1] x [0, m - 1] asserts
 * q updates create O(q log n log m) nodes from a pool of MAXNODES, 16 bytes each for long long, about 205 MB by default
 * The pool is reserved up front, but only the pages of created nodes are touched; running out asserts, raise MAXNODES
 * The row roots take 4 * (n + 1) bytes, so n is bounded by memory and m only by the pool
 * fenwick_tree_2D_sparse.cpp does point add and rectangle sum on a fixed hash table, 268 MB by default, this one pays only for created nodes
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
    /// Two overlapping rectangle adds
    FenwickImplicit2D<long long> fen(1000, 1000);
    fen.update(2, 4, 5, 6, 10);
    fen.update(3, 5, 7, 8, 5);
    assert(fen.query(2, 4) == 10);
    assert(fen.query(4, 5) == 15);
    assert(fen.query(6, 7) == 5);
    assert(fen.query(1, 3) == 0);

    /// 2e9 columns, only the created nodes take memory
    FenwickImplicit2D<long long> wide(100000, 2000000000);
    wide.update(4, 6, 99999, 1999999999, 3);
    assert(wide.query(99999, 1999999999) == 3);
    assert(wide.query(4, 5) == 0);

    return 0;
}
