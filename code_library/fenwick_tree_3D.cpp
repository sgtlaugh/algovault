/***
 *
 * 3D Fenwick Tree (Binary Indexed Tree)
 * Supports point update + range query, range update + point query, and range update + range query variants
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
 * Point updates, range queries
 * update(x, y, z, v): Add v to index [x, y, z]
 * query(x1, y1, z1, x2, y2, z2): Return sum of cube from [x1,y1,z1] to [x2,y2,z2]
***/

template <typename T>
struct FenwickPointUpdate3D{
    int n, m, r;
    vector<vector<vector<T>>> tree;

    FenwickPointUpdate3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0))) {}

    void update(int x, int y, int z, T v){
        for (int i = x; i <= n; i += i & -i){
            for (int j = y; j <= m; j += j & -j){
                for (int k = z; k <= r; k += k & -k){
                    tree[i][j][k] += v;
                }
            }
        }
    }

    T query(int x, int y, int z){
        if (x <= 0 || y <= 0 || z <= 0) return 0;
        T res = 0;
        for (int i = x; i > 0; i -= i & -i){
            for (int j = y; j > 0; j -= j & -j){
                for (int k = z; k > 0; k -= k & -k){
                    res += tree[i][j][k];
                }
            }
        }
        return res;
    }

    T query(int x1, int y1, int z1, int x2, int y2, int z2){
        if (x1 > x2 || y1 > y2 || z1 > z2) return 0;
        return query(x2, y2, z2)
             - query(x1 - 1, y2, z2) - query(x2, y1 - 1, z2) - query(x2, y2, z1 - 1)
             + query(x1 - 1, y1 - 1, z2) + query(x1 - 1, y2, z1 - 1) + query(x2, y1 - 1, z1 - 1)
             - query(x1 - 1, y1 - 1, z1 - 1);
    }
};

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

/***
 * Range updates, range queries
 * update(x1, y1, z1, x2, y2, z2, v): Add v to cube from [x1,y1,z1] to [x2,y2,z2]
 * query(x1, y1, z1, x2, y2, z2): Return sum of cube from [x1,y1,z1] to [x2,y2,z2]
***/

template <typename T>
struct FenwickFull3D{
    int n, m, r;
    vector<vector<vector<T>>> tree[8];

    FenwickFull3D(int n = 0, int m = 0, int r = 0) : n(n), m(m), r(r) {
        for (int i = 0; i < 8; i++){
            tree[i].assign(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0)));
        }
    }

    void update_base(int p, int q, int s, T v){
        if (p <= 0 || q <= 0 || s <= 0 || p > n || q > m || s > r) return;

        T a = p - 1, b = q - 1, c = s - 1;
        for (int i = p; i <= n; i += i & -i){
            for (int j = q; j <= m; j += j & -j){
                for (int k = s; k <= r; k += k & -k){
                    tree[0][i][j][k] += v;
                    tree[1][i][j][k] += v * a;
                    tree[2][i][j][k] += v * b;
                    tree[3][i][j][k] += v * c;
                    tree[4][i][j][k] += v * a * b;
                    tree[5][i][j][k] += v * a * c;
                    tree[6][i][j][k] += v * b * c;
                    tree[7][i][j][k] += v * a * b * c;
                }
            }
        }
    }

    void update(int x1, int y1, int z1, int x2, int y2, int z2, T v){
        if (x1 > x2 || y1 > y2 || z1 > z2) return;

        update_base(x1, y1, z1, v);
        update_base(x2 + 1, y1, z1, -v);
        update_base(x1, y2 + 1, z1, -v);
        update_base(x1, y1, z2 + 1, -v);
        update_base(x2 + 1, y2 + 1, z1, v);
        update_base(x2 + 1, y1, z2 + 1, v);
        update_base(x1, y2 + 1, z2 + 1, v);
        update_base(x2 + 1, y2 + 1, z2 + 1, -v);
    }

    T query_base(int p, int q, int s){
        if (p <= 0 || q <= 0 || s <= 0) return 0;

        T t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
        for (int i = p; i > 0; i -= i & -i){
            for (int j = q; j > 0; j -= j & -j){
                for (int k = s; k > 0; k -= k & -k){
                    t0 += tree[0][i][j][k];
                    t1 += tree[1][i][j][k];
                    t2 += tree[2][i][j][k];
                    t3 += tree[3][i][j][k];
                    t4 += tree[4][i][j][k];
                    t5 += tree[5][i][j][k];
                    t6 += tree[6][i][j][k];
                    t7 += tree[7][i][j][k];
                }
            }
        }
        return t0 * p * q * s - t1 * q * s - t2 * p * s - t3 * p * q
             + t4 * s + t5 * q + t6 * p - t7;
    }

    T query(int x1, int y1, int z1, int x2, int y2, int z2){
        if (x1 > x2 || y1 > y2 || z1 > z2) return 0;
        return query_base(x2, y2, z2)
             - query_base(x1 - 1, y2, z2) - query_base(x2, y1 - 1, z2) - query_base(x2, y2, z1 - 1)
             + query_base(x1 - 1, y1 - 1, z2) + query_base(x1 - 1, y2, z1 - 1) + query_base(x2, y1 - 1, z1 - 1)
             - query_base(x1 - 1, y1 - 1, z1 - 1);
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate3D<int>(10, 10, 10);
    fen.update(2, 3, 4, 5);
    fen.update(7, 8, 9, 10);
    fen.update(5, 6, 7, 3);

    assert(fen.query(2, 3, 4, 2, 3, 4) == 5);
    assert(fen.query(1, 1, 1, 10, 10, 10) == 18);
    assert(fen.query(5, 6, 7, 8, 9, 10) == 13);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate3D<int>(10, 10, 10);
    fen2.update(2, 3, 4, 5, 6, 7, 1);
    fen2.update(4, 5, 6, 8, 9, 10, 1);

    assert(fen2.query(3, 4, 5) == 1);
    assert(fen2.query(5, 6, 7) == 2);
    assert(fen2.query(7, 8, 9) == 1);
    assert(fen2.query(1, 1, 1) == 0);

    /// Range update, range query
    auto fen3 = FenwickFull3D<long long>(10, 10, 10);
    fen3.update(2, 3, 4, 5, 6, 7, 2);
    fen3.update(4, 5, 6, 8, 9, 10, 3);

    assert(fen3.query(2, 3, 4, 5, 6, 7) == 152);
    assert(fen3.query(4, 5, 6, 8, 9, 10) == 391);
    assert(fen3.query(4, 5, 6, 5, 6, 7) == 40);

    return 0;
}
