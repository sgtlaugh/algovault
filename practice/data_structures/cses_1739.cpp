/***
 *
 * 2D Fenwick Tree (Binary Indexed Tree)
 * Sums over a grid under updates: point add and rectangle sum, rectangle add and cell value, or rectangle add and rectangle sum
 *
 * Complexity:
 *   - build: O(n m)
 *   - update, prefix, query: O(log n log m)
 *
 * FenwickPointUpdate2D<T> fen(n, m): point add, rectangle sum, over an n x m grid of zeros
 * fen.update(i, j, v): add v to a[i][j]
 * fen.prefix(i, j): sum of a[0..i][0..j], 0 when i or j is -1
 * fen.query(i, j, k, l): sum of a[i..k][j..l], 0 when i > k or j > l
 *
 * FenwickRangeUpdate2D<T> fen(n, m): rectangle add, cell value, over an n x m grid of zeros
 * fen.update(i, j, k, l, v): add v to a[i..k][j..l], nothing when i > k or j > l
 * fen.query(i, j): a[i][j]
 *
 * FenwickFull2D<T> fen(n, m): rectangle add, rectangle sum, over an n x m grid of zeros
 * fen.update(i, j, k, l, v): add v to a[i..k][j..l], nothing when i > k or j > l
 * fen.prefix(i, j): sum of a[0..i][0..j], 0 when i or j is -1
 * fen.query(i, j, k, l): sum of a[i..k][j..l], 0 when i > k or j > l
 *
 * Cells are 0-based, a cell outside [0, n - 1] x [0, m - 1] asserts
 * FenwickFull2D keeps four values per cell, four times the memory of the other two
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickPointUpdate2D{
    int n, m;
    vector<vector<T>> tree;

    FenwickPointUpdate2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<T>(m + 1, 0)) {}

    void update(int i, int j, T v){
        assert(0 <= i && i < n && 0 <= j && j < m);
        for (int x = i + 1; x <= n; x += x & -x){
            for (int y = j + 1; y <= m; y += y & -y){
                tree[x][y] += v;
            }
        }
    }

    T prefix(int i, int j){
        assert(-1 <= i && i < n && -1 <= j && j < m);
        T res = 0;
        for (int x = i + 1; x > 0; x -= x & -x){
            for (int y = j + 1; y > 0; y -= y & -y){
                res += tree[x][y];
            }
        }
        return res;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        assert(0 <= i && k < n && 0 <= j && l < m);
        return prefix(k, l) - prefix(i - 1, l) - prefix(k, j - 1) + prefix(i - 1, j - 1);
    }
};

template <typename T>
struct FenwickRangeUpdate2D{
    int n, m;
    vector<vector<T>> tree;

    FenwickRangeUpdate2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<T>(m + 1, 0)) {}

    void update(int i, int j, int k, int l, T v){
        if (i > k || j > l) return;
        assert(0 <= i && k < n && 0 <= j && l < m);
        add_suffix(i, j, v);
        add_suffix(k + 1, j, -v);
        add_suffix(k + 1, l + 1, v);
        add_suffix(i, l + 1, -v);
    }

    T query(int i, int j){
        assert(0 <= i && i < n && 0 <= j && j < m);
        T res = 0;
        for (int x = i + 1; x > 0; x -= x & -x){
            for (int y = j + 1; y > 0; y -= y & -y){
                res += tree[x][y];
            }
        }
        return res;
    }

    void add_suffix(int i, int j, T v){
        for (int x = i + 1; x <= n; x += x & -x){
            for (int y = j + 1; y <= m; y += y & -y){
                tree[x][y] += v;
            }
        }
    }
};

template <typename T>
struct FenwickFull2D{
    int n, m;
    vector<vector<array<T, 4>>> tree;  /// Sums of v, v*b, v*a, v*a*b per node; one array shares a cache line, about 2x faster than four trees

    FenwickFull2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<array<T, 4>>(m + 1, array<T, 4>{})) {}

    void update(int i, int j, int k, int l, T v){
        if (i > k || j > l) return;
        assert(0 <= i && k < n && 0 <= j && l < m);
        add_suffix(i, j, v);
        add_suffix(k + 1, j, -v);
        add_suffix(k + 1, l + 1, v);
        add_suffix(i, l + 1, -v);
    }

    /// A difference v at (a, b) adds v * (p - a) * (q - b) = v*p*q - v*b*p - v*a*q + v*a*b to the sum of the top-left p x q cells
    T prefix(int p, int q){
        assert(-1 <= p && p < n && -1 <= q && q < m);
        p++, q++;
        T v = 0, vb = 0, va = 0, vab = 0;
        for (int i = p; i > 0; i -= i & -i){
            for (int j = q; j > 0; j -= j & -j){
                auto& t = tree[i][j];
                v += t[0];
                vb += t[1];
                va += t[2];
                vab += t[3];
            }
        }
        return v * p * q - vb * p - va * q + vab;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        assert(0 <= i && k < n && 0 <= j && l < m);
        return prefix(k, l) - prefix(i - 1, l) - prefix(k, j - 1) + prefix(i - 1, j - 1);
    }

    /// The products stay inside the loop: hoisting them spills registers and measured 5% slower
    void add_suffix(int a, int b, T v){
        for (int i = a + 1; i <= n; i += i & -i){
            for (int j = b + 1; j <= m; j += j & -j){
                auto& t = tree[i][j];
                t[0] += v;
                t[1] += v * b;
                t[2] += v * a;
                t[3] += v * a * b;
            }
        }
    }
};

int main(){
    int n, q;
    scanf("%d %d", &n, &q);

    FenwickPointUpdate2D<int> fen(n, n);
    vector<string> grid(n);
    char buf[1005];
    for (int y = 0; y < n; y++){
        scanf("%s", buf);
        grid[y] = buf;
        for (int x = 0; x < n; x++){
            if (grid[y][x] == '*') fen.update(y, x, 1);
        }
    }

    while (q--){
        int t;
        scanf("%d", &t);
        if (t == 1){
            int y, x;
            scanf("%d %d", &y, &x);
            y--, x--;
            fen.update(y, x, grid[y][x] == '*' ? -1 : 1);
            grid[y][x] = grid[y][x] == '*' ? '.' : '*';
        }
        else{
            int y1, x1, y2, x2;
            scanf("%d %d %d %d", &y1, &x1, &y2, &x2);
            printf("%d\n", fen.query(y1 - 1, x1 - 1, y2 - 1, x2 - 1));
        }
    }
    return 0;
}
