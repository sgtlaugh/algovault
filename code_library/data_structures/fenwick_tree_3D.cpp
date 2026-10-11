/***
 *
 * 3D Fenwick Tree (Binary Indexed Tree)
 * Supports point update + range query, range update + point query, and range update + range query variants
 * 0-based cells, boxes [x1, x2] x [y1, y2] x [z1, z2] inclusive, a cell outside the n x m x r grid asserts
 * A box with x1 > x2, y1 > y2 or z1 > z2 is empty: updates on it do nothing and query returns 0
 *
 * Complexity:
 *   - O(n * m * r) to build
 *   - O(log n * log m * log r) per update/query
 *
 * FenwickFull3D<long long> f(n, m, r);
 * f.update(x1, y1, z1, x2, y2, z2, v); f.query(x1, y1, z1, x2, y2, z2);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Point updates, range queries
 * update(x, y, z, v): Add v to cell [x, y, z]
 * prefix(x, y, z): Return sum of cube from [0,0,0] to [x,y,z], 0 when any coordinate is -1
 * query(x1, y1, z1, x2, y2, z2): Return sum of cube from [x1,y1,z1] to [x2,y2,z2]
***/

template <typename T>
struct FenwickPointUpdate3D{
    int n, m, r;
    vector<vector<vector<T>>> tree;

    FenwickPointUpdate3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0))) {}

    void update(int x, int y, int z, T v){
        assert(0 <= x && x < n && 0 <= y && y < m && 0 <= z && z < r);
        for (int i = x + 1; i <= n; i += i & -i){
            for (int j = y + 1; j <= m; j += j & -j){
                for (int k = z + 1; k <= r; k += k & -k){
                    tree[i][j][k] += v;
                }
            }
        }
    }

    T prefix(int x, int y, int z){
        assert(-1 <= x && x < n && -1 <= y && y < m && -1 <= z && z < r);
        T res = 0;
        for (int i = x + 1; i > 0; i -= i & -i){
            for (int j = y + 1; j > 0; j -= j & -j){
                for (int k = z + 1; k > 0; k -= k & -k){
                    res += tree[i][j][k];
                }
            }
        }
        return res;
    }

    T query(int x1, int y1, int z1, int x2, int y2, int z2){
        if (x1 > x2 || y1 > y2 || z1 > z2) return 0;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m && 0 <= z1 && z2 < r);
        return prefix(x2, y2, z2)
             - prefix(x1 - 1, y2, z2) - prefix(x2, y1 - 1, z2) - prefix(x2, y2, z1 - 1)
             + prefix(x1 - 1, y1 - 1, z2) + prefix(x1 - 1, y2, z1 - 1) + prefix(x2, y1 - 1, z1 - 1)
             - prefix(x1 - 1, y1 - 1, z1 - 1);
    }
};

/***
 * Range updates, point queries
 * update(x1, y1, z1, x2, y2, z2, v): Add v to cube from [x1,y1,z1] to [x2,y2,z2]
 * query(x, y, z): Return value at cell [x, y, z]
***/

template <typename T>
struct FenwickRangeUpdate3D{
    int n, m, r;
    vector<vector<vector<T>>> tree;

    FenwickRangeUpdate3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0))) {}

    void update(int x1, int y1, int z1, int x2, int y2, int z2, T v){
        if (x1 > x2 || y1 > y2 || z1 > z2) return;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m && 0 <= z1 && z2 < r);

        add_prefix(x2, y2, z2, v);
        add_prefix(x1 - 1, y1 - 1, z2, v);
        add_prefix(x1 - 1, y2, z1 - 1, v);
        add_prefix(x2, y1 - 1, z1 - 1, v);
        add_prefix(x1 - 1, y2, z2, -v);
        add_prefix(x2, y1 - 1, z2, -v);
        add_prefix(x2, y2, z1 - 1, -v);
        add_prefix(x1 - 1, y1 - 1, z1 - 1, -v);
    }

    T query(int x, int y, int z){
        assert(0 <= x && x < n && 0 <= y && y < m && 0 <= z && z < r);
        T res = 0;
        for (int i = x + 1; i <= n; i += i & -i){
            for (int j = y + 1; j <= m; j += j & -j){
                for (int k = z + 1; k <= r; k += k & -k){
                    res += tree[i][j][k];
                }
            }
        }
        return res;
    }

private:
    void add_prefix(int x, int y, int z, T v){
        for (int i = x + 1; i > 0; i -= i & -i){
            for (int j = y + 1; j > 0; j -= j & -j){
                for (int k = z + 1; k > 0; k -= k & -k){
                    tree[i][j][k] += v;
                }
            }
        }
    }
};

/***
 * Range updates, range queries
 * update(x1, y1, z1, x2, y2, z2, v): Add v to cube from [x1,y1,z1] to [x2,y2,z2]
 * prefix(x, y, z): Return sum of cube from [0,0,0] to [x,y,z], 0 when any coordinate is -1
 * query(x1, y1, z1, x2, y2, z2): Return sum of cube from [x1,y1,z1] to [x2,y2,z2]
***/

template <typename T>
struct FenwickFull3D{
    int n, m, r;
    vector<vector<vector<array<T, 8>>>> tree;  /// A cell's eight values stay adjacent, about 2x faster than eight trees

    FenwickFull3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<array<T, 8>>>(m + 1, vector<array<T, 8>>(r + 1, array<T, 8>{}))) {}

    void update(int x1, int y1, int z1, int x2, int y2, int z2, T v){
        if (x1 > x2 || y1 > y2 || z1 > z2) return;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m && 0 <= z1 && z2 < r);

        update_base(x1, y1, z1, v);
        update_base(x2 + 1, y1, z1, -v);
        update_base(x1, y2 + 1, z1, -v);
        update_base(x1, y1, z2 + 1, -v);
        update_base(x2 + 1, y2 + 1, z1, v);
        update_base(x2 + 1, y1, z2 + 1, v);
        update_base(x1, y2 + 1, z2 + 1, v);
        update_base(x2 + 1, y2 + 1, z2 + 1, -v);
    }

    T prefix(int p, int q, int s){
        assert(-1 <= p && p < n && -1 <= q && q < m && -1 <= s && s < r);
        p++, q++, s++;
        T t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
        for (int i = p; i > 0; i -= i & -i){
            for (int j = q; j > 0; j -= j & -j){
                for (int k = s; k > 0; k -= k & -k){
                    auto& t = tree[i][j][k];
                    t0 += t[0];
                    t1 += t[1];
                    t2 += t[2];
                    t3 += t[3];
                    t4 += t[4];
                    t5 += t[5];
                    t6 += t[6];
                    t7 += t[7];
                }
            }
        }
        return t0 * p * q * s - t1 * q * s - t2 * p * s - t3 * p * q
             + t4 * s + t5 * q + t6 * p - t7;
    }

    T query(int x1, int y1, int z1, int x2, int y2, int z2){
        if (x1 > x2 || y1 > y2 || z1 > z2) return 0;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m && 0 <= z1 && z2 < r);
        return prefix(x2, y2, z2)
             - prefix(x1 - 1, y2, z2) - prefix(x2, y1 - 1, z2) - prefix(x2, y2, z1 - 1)
             + prefix(x1 - 1, y1 - 1, z2) + prefix(x1 - 1, y2, z1 - 1) + prefix(x2, y1 - 1, z1 - 1)
             - prefix(x1 - 1, y1 - 1, z1 - 1);
    }

private:
    void update_base(int p, int q, int s, T v){
        T a = p, b = q, c = s;
        for (int i = p + 1; i <= n; i += i & -i){
            for (int j = q + 1; j <= m; j += j & -j){
                for (int k = s + 1; k <= r; k += k & -k){
                    auto& t = tree[i][j][k];
                    t[0] += v;
                    t[1] += v * a;
                    t[2] += v * b;
                    t[3] += v * c;
                    t[4] += v * a * b;
                    t[5] += v * a * c;
                    t[6] += v * b * c;
                    t[7] += v * a * b * c;
                }
            }
        }
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate3D<int>(10, 10, 10);
    fen.update(1, 2, 3, 5);
    fen.update(6, 7, 8, 10);
    fen.update(4, 5, 6, 3);

    assert(fen.query(1, 2, 3, 1, 2, 3) == 5);
    assert(fen.query(0, 0, 0, 9, 9, 9) == 18);
    assert(fen.query(4, 5, 6, 7, 8, 9) == 13);
    assert(fen.prefix(4, 5, 6) == 8);
    assert(fen.prefix(9, -1, 9) == 0);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate3D<int>(10, 10, 10);
    fen2.update(1, 2, 3, 4, 5, 6, 1);
    fen2.update(3, 4, 5, 7, 8, 9, 1);

    assert(fen2.query(2, 3, 4) == 1);
    assert(fen2.query(4, 5, 6) == 2);
    assert(fen2.query(6, 7, 8) == 1);
    assert(fen2.query(7, 8, 9) == 1);
    assert(fen2.query(0, 0, 0) == 0);

    /// Range update, range query
    auto fen3 = FenwickFull3D<long long>(10, 10, 10);
    fen3.update(1, 2, 3, 4, 5, 6, 2);
    fen3.update(3, 4, 5, 7, 8, 9, 3);

    assert(fen3.query(1, 2, 3, 4, 5, 6) == 152);
    assert(fen3.query(3, 4, 5, 7, 8, 9) == 391);
    assert(fen3.query(3, 4, 5, 4, 5, 6) == 40);
    assert(fen3.prefix(2, 3, 4) == 16);
    assert(fen3.prefix(-1, 9, 9) == 0);

    return 0;
}
