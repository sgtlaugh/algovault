/***
 *
 * 2D Fenwick Tree with Implicit Segment Trees for large matrices with sparse updates
 * Each Fenwick row uses implicit segment tree for columns - allocates nodes on demand
 * Range updates and point queries only
 *
 * Complexity:
 *   - O(n) to construct, for the row roots
 *   - O(log^2 n) per update/query, a range update is 4 point updates of O(log n) rows each creating up to log n + 1 nodes
 *   - O(q log^2 n) nodes for q updates
 *
 * Nodes come from a pool of MAXNODES (16 bytes each for long long, ~205 MB by default), each range update takes up to 4 * log^2(N)
 * The pool is reserved up front but only memory for nodes actually created is touched
 * The row roots take another 4 * (N + 1) bytes, so N is bounded by memory
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
    static const int MAXNODES = 12800000; // Adjust based on memory constraints

    struct Node{
        T sum;
        int l, r;
    };

    int n;
    vector<int> root;
    vector<Node> nodes;

    /// reserve without filling, so pool pages are only touched once nodes are used, node 0 is the empty child
    FenwickImplicit2D(int n = 0) : n(n), root(n + 1, 0){
        nodes.reserve(MAXNODES);
        nodes.push_back({0, 0, 0});
    }

    int update_seg(int cur, int a, int b, int p, T v){
        if (!cur){
            assert((int)nodes.size() < MAXNODES);  /// out of pool nodes, raise MAXNODES
            cur = nodes.size();
            nodes.push_back({0, 0, 0});
        }
        nodes[cur].sum += v;

        if (a != b){
            int m = (a + b) >> 1;
            if (p <= m) nodes[cur].l = update_seg(nodes[cur].l, a, m, p, v);
            else nodes[cur].r = update_seg(nodes[cur].r, m + 1, b, p, v);
        }
        return cur;
    }

    T query_seg(int cur, int a, int b, int r){
        if (!cur) return 0;
        if (r >= b) return nodes[cur].sum;

        int m = (a + b) >> 1;
        if (r <= m) return query_seg(nodes[cur].l, a, m, r);
        return nodes[nodes[cur].l].sum + query_seg(nodes[cur].r, m + 1, b, r);
    }

    void update_fen(int x, int y, T v){
        if (y > n) return;  /// column n + 1 from a range update ending at n, otherwise it lands on leaf n
        for (int i = x; i <= n; i += i & -i){
            root[i] = update_seg(root[i], 1, n, y, v);
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
