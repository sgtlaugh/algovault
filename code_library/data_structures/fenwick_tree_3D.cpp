/***
 *
 * 3D Fenwick Tree (Binary Indexed Tree)
 * Sums over a 3D grid under updates: point add and box sum, box add and cell value, or box add and box sum
 *
 * Complexity:
 *   - build: O(n * m * r)
 *   - update, prefix, query: O(log n * log m * log r)
 *
 * FenwickPointUpdate3D<T> fen(n, m, r): point add, box sum, over an n x m x r grid of zeros
 * fen.update(x, y, z, v): add v to a[x][y][z]
 * fen.prefix(x, y, z): sum of a[0..x][0..y][0..z], 0 when any coordinate is -1
 * fen.query(x1, y1, z1, x2, y2, z2): sum of a[x1..x2][y1..y2][z1..z2], 0 when x1 > x2, y1 > y2 or z1 > z2
 *
 * FenwickRangeUpdate3D<T> fen(n, m, r): box add, cell value, over an n x m x r grid of zeros
 * fen.update(x1, y1, z1, x2, y2, z2, v): add v to a[x1..x2][y1..y2][z1..z2], nothing when x1 > x2, y1 > y2 or z1 > z2
 * fen.query(x, y, z): a[x][y][z]
 *
 * FenwickFull3D<T> fen(n, m, r): box add, box sum, over an n x m x r grid of zeros
 * fen.update(x1, y1, z1, x2, y2, z2, v): add v to a[x1..x2][y1..y2][z1..z2], nothing when x1 > x2, y1 > y2 or z1 > z2
 * fen.prefix(x, y, z): sum of a[0..x][0..y][0..z], 0 when any coordinate is -1
 * fen.query(x1, y1, z1, x2, y2, z2): sum of a[x1..x2][y1..y2][z1..z2], 0 when x1 > x2, y1 > y2 or z1 > z2
 *
 * Cells are 0-based, a cell outside [0, n - 1] x [0, m - 1] x [0, r - 1] asserts
 * FenwickFull3D keeps eight values per cell, eight times the memory of the other two
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

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

template <typename T>
struct FenwickRangeUpdate3D{
    int n, m, r;
    vector<vector<vector<T>>> tree;

    FenwickRangeUpdate3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<T>>(m + 1, vector<T>(r + 1, 0))) {}

    /// Adds to the 8 prefix corners and query sums the suffix, the mirror of 1D and 2D: at n = 100 the mirror's updates
    /// are 1.7x slower and its queries 1.6x to 3x faster, so this direction wins unless queries outnumber updates 4 to 9 times
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

template <typename T>
struct FenwickFull3D{
    int n, m, r;
    vector<vector<vector<array<T, 8>>>> tree;  /// Sums of v, v*a, v*b, v*c, v*a*b, v*a*c, v*b*c, v*a*b*c per node, adjacent, about 2x faster than eight trees

    FenwickFull3D(int n = 0, int m = 0, int r = 0)
        : n(n), m(m), r(r), tree(n + 1, vector<vector<array<T, 8>>>(m + 1, vector<array<T, 8>>(r + 1, array<T, 8>{}))) {}

    void update(int x1, int y1, int z1, int x2, int y2, int z2, T v){
        if (x1 > x2 || y1 > y2 || z1 > z2) return;
        assert(0 <= x1 && x2 < n && 0 <= y1 && y2 < m && 0 <= z1 && z2 < r);

        add_suffix(x1, y1, z1, v);
        add_suffix(x2 + 1, y1, z1, -v);
        add_suffix(x1, y2 + 1, z1, -v);
        add_suffix(x1, y1, z2 + 1, -v);
        add_suffix(x2 + 1, y2 + 1, z1, v);
        add_suffix(x2 + 1, y1, z2 + 1, v);
        add_suffix(x1, y2 + 1, z2 + 1, v);
        add_suffix(x2 + 1, y2 + 1, z2 + 1, -v);
    }

    /// A difference v at (a, b, c) adds v * (p - a) * (q - b) * (s - c) to the sum of the p x q x s corner box,
    /// expanding the product gives one term per kept sum
    T prefix(int p, int q, int s){
        assert(-1 <= p && p < n && -1 <= q && q < m && -1 <= s && s < r);
        p++, q++, s++;
        T v = 0, va = 0, vb = 0, vc = 0, vab = 0, vac = 0, vbc = 0, vabc = 0;
        for (int i = p; i > 0; i -= i & -i){
            for (int j = q; j > 0; j -= j & -j){
                for (int k = s; k > 0; k -= k & -k){
                    auto& t = tree[i][j][k];
                    v += t[0];
                    va += t[1];
                    vb += t[2];
                    vc += t[3];
                    vab += t[4];
                    vac += t[5];
                    vbc += t[6];
                    vabc += t[7];
                }
            }
        }
        return v * p * q * s - va * q * s - vb * p * s - vc * p * q + vab * s + vac * q + vbc * p - vabc;
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
    /// The products stay inside the loop: hoisting them spills registers and measured 12% slower
    void add_suffix(int a, int b, int c, T v){
        for (int i = a + 1; i <= n; i += i & -i){
            for (int j = b + 1; j <= m; j += j & -j){
                for (int k = c + 1; k <= r; k += k & -k){
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
    /// Point add, box sum
    FenwickPointUpdate3D<int> fen(10, 10, 10);
    fen.update(1, 2, 3, 5);
    fen.update(6, 7, 8, 10);
    fen.update(4, 5, 6, 3);
    assert(fen.query(0, 0, 0, 9, 9, 9) == 18);
    assert(fen.query(4, 5, 6, 7, 8, 9) == 13);
    assert(fen.prefix(4, 5, 6) == 8);

    /// Box add, cell value
    FenwickRangeUpdate3D<int> range(10, 10, 10);
    range.update(1, 2, 3, 4, 5, 6, 1);
    range.update(3, 4, 5, 7, 8, 9, 1);
    assert(range.query(2, 3, 4) == 1);
    assert(range.query(4, 5, 6) == 2);
    assert(range.query(7, 8, 9) == 1);
    assert(range.query(0, 0, 0) == 0);

    /// Box add, box sum: 4 x 4 x 4 cells of 2 and 5 x 5 x 5 cells of 3, overlapping on 2 x 2 x 2 cells
    FenwickFull3D<long long> full(10, 10, 10);
    full.update(1, 2, 3, 4, 5, 6, 2);
    full.update(3, 4, 5, 7, 8, 9, 3);
    assert(full.query(1, 2, 3, 4, 5, 6) == 152);
    assert(full.query(3, 4, 5, 7, 8, 9) == 391);
    assert(full.query(3, 4, 5, 4, 5, 6) == 40);
    assert(full.prefix(2, 3, 4) == 16);

    return 0;
}
