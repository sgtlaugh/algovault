/***
 *
 * 2D Fenwick Tree (Binary Indexed Tree)
 * Supports point update + range query, range update + point query, and range update + range query variants
 * 0-based cells, rectangles [i, k] x [j, l] inclusive, a cell outside [0, n - 1] x [0, m - 1] asserts
 * i > k or j > l is an empty rectangle: updates on it do nothing and query returns 0
 *
 * Complexity:
 *   - O(n * m) to build
 *   - O(log n * log m) per update/query
 *
 * FenwickFull2D<long long> f(n, m);
 * f.update(i, j, k, l, v); f.query(i, j, k, l);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Point updates, range queries
 * update(i, j, v): Add v to cell [i, j]
 * prefix(i, j): Return sum of rectangle from [0, 0] to [i, j], 0 when i or j is -1
 * query(i, j, k, l): Return sum of rectangle from [i, j] to [k, l]
***/

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

/***
 * Range updates, point queries
 * update(i, j, k, l, v): Add v to rectangle from [i, j] to [k, l]
 * query(i, j): Return value at cell [i, j]
***/

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

private:
    void add_suffix(int i, int j, T v){
        for (int x = i + 1; x <= n; x += x & -x){
            for (int y = j + 1; y <= m; y += y & -y){
                tree[x][y] += v;
            }
        }
    }
};

/***
 * Range updates, range queries
 * update(i, j, k, l, v): Add v to rectangle from [i, j] to [k, l]
 * prefix(i, j): Return sum of rectangle from [0, 0] to [i, j], 0 when i or j is -1
 * query(i, j, k, l): Return sum of rectangle from [i, j] to [k, l]
***/

template <typename T>
struct FenwickFull2D{
    int n, m;
    vector<vector<array<T, 4>>> tree;  /// A cell's four values share a cache line, about 2x faster than four trees

    FenwickFull2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<array<T, 4>>(m + 1, array<T, 4>{})) {}

    void update(int i, int j, int k, int l, T v){
        if (i > k || j > l) return;
        assert(0 <= i && k < n && 0 <= j && l < m);
        add_suffix(i, j, v);
        add_suffix(k + 1, j, -v);
        add_suffix(k + 1, l + 1, v);
        add_suffix(i, l + 1, -v);
    }

    T prefix(int p, int q){
        assert(-1 <= p && p < n && -1 <= q && q < m);
        p++, q++;
        T x = 0, y = 0, z = 0;
        for (int i = p; i > 0; i -= i & -i){
            T c = 0, d = 0;
            for (int j = q; j > 0; j -= j & -j){
                auto& t = tree[i][j];
                c += t[0];
                d += t[1];
                y += t[2];
                z += t[3];
            }
            x += c * q - d;
        }
        return x * p - y * q + z;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        assert(0 <= i && k < n && 0 <= j && l < m);
        return prefix(k, l) - prefix(i - 1, l) - prefix(k, j - 1) + prefix(i - 1, j - 1);
    }

private:
    void add_suffix(int p, int q, T v){
        T c = p, d = q;
        for (int i = p + 1; i <= n; i += i & -i){
            for (int j = q + 1; j <= m; j += j & -j){
                auto& t = tree[i][j];
                t[0] += v;
                t[1] += v * d;
                t[2] += v * c;
                t[3] += v * c * d;
            }
        }
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate2D<int>(10, 10);
    fen.update(1, 2, 5);
    fen.update(4, 6, 3);
    fen.update(3, 3, 7);
    fen.update(9, 9, 1);

    assert(fen.query(1, 2, 1, 2) == 5);
    assert(fen.query(0, 0, 9, 9) == 16);
    assert(fen.query(3, 3, 5, 7) == 10);
    assert(fen.query(4, 4, 3, 3) == 0);
    assert(fen.prefix(3, 3) == 12);
    assert(fen.prefix(-1, 9) == 0);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate2D<int>(10, 10);
    fen2.update(1, 2, 4, 6, 4);
    fen2.update(3, 3, 9, 9, 2);

    assert(fen2.query(1, 2) == 4);
    assert(fen2.query(4, 5) == 6);
    assert(fen2.query(3, 3) == 6);
    assert(fen2.query(0, 0) == 0);
    assert(fen2.query(9, 9) == 2);

    /// Range update, range query
    auto fen3 = FenwickFull2D<long long>(10, 10);
    fen3.update(1, 2, 4, 6, 4);
    fen3.update(3, 3, 9, 9, 2);

    assert(fen3.query(1, 2, 4, 6) == 96);
    assert(fen3.query(3, 3, 9, 9) == 130);
    assert(fen3.query(1, 2, 9, 9) == 178);
    assert(fen3.prefix(1, 2) == 4);

    return 0;
}
