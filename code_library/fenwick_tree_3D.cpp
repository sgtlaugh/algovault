/***
 *
 * 3D Fenwick Tree (Binary Indexed Tree)
 * Range updates and point queries
 * All operations in O(log n * log m * log r)
 * 1-based indexing for elements
 *
 * Complexity:
 *   - O(n * m * r) to build
 *   - O(log n * log m * log r) per update/query
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Range updates, point queries
 * update(x1, y1, z1, x2, y2, z2, v): Add v to cube from [x1,y1,z1] to [x2,y2,z2]
 * query(x, y, z): Return value at index [x, y, z]
***/

template <typename T>
struct FenwickRangeUpdate3D{
    int n, m, r;
    vector<vector<vector<T>>> tree;

    FenwickRangeUpdate3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0))) {}

    void update(int x, int y, int z, T v){
        if (x < 0 || y < 0 || z < 0 || x > n || y > m || z > r) return;

        for (int i = x; i > 0; i -= i & -i){
            for (int j = y; j > 0; j -= j & -j){
                for (int k = z; k > 0; k -= k & -k){
                    tree[i][j][k] += v;
                }
            }
        }
    }

    void update(int x1, int y1, int z1, int x2, int y2, int z2, T v){
        if (x1 > x2 || y1 > y2 || z1 > z2) return;

        update(x2, y2, z2, v);
        update(x1 - 1, y1 - 1, z2, v);
        update(x1 - 1, y2, z1 - 1, v);
        update(x2, y1 - 1, z1 - 1, v);
        update(x1 - 1, y2, z2, -v);
        update(x2, y1 - 1, z2, -v);
        update(x2, y2, z1 - 1, -v);
        update(x1 - 1, y1 - 1, z1 - 1, -v);
    }

    T query(int x, int y, int z){
        T res = 0;
        for (int i = x; i <= n; i += i & -i){
            for (int j = y; j <= m; j += j & -j){
                for (int k = z; k <= r; k += k & -k){
                    res += tree[i][j][k];
                }
            }
        }
        return res;
    }
};

int main(){
    /// Test 3D range updates
    auto fen = FenwickRangeUpdate3D<int>(10, 10, 10);

    fen.update(2, 3, 4, 5, 6, 7, 1);
    assert(fen.query(3, 4, 5) == 1);
    assert(fen.query(6, 7, 8) == 0);
    assert(fen.query(1, 1, 1) == 0);

    fen.update(4, 5, 6, 8, 9, 10, 1);
    assert(fen.query(5, 6, 7) == 2);  // Overlap: both updates
    assert(fen.query(7, 8, 9) == 1);  // Second update only
    assert(fen.query(3, 4, 5) == 1);  // First update only

    return 0;
}
