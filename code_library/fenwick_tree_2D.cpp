/***
 *
 * 2D Fenwick Tree (Binary Indexed Tree)
 * Supports point update + range query, range update + point query, and range update + range query variants
 * All operations in O(log n * log m)
 * 1-based indexing for elements
 *
 * Complexity:
 *   - O(n * m) to build
 *   - O(log n * log m) per update/query
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Point updates, range queries
 * update(i, j, v): Add v to index [i, j]
 * query(i, j, k, l): Return sum of rectangle from [i, j] to [k, l]
***/

template <typename T>
struct FenwickPointUpdate2D{
    int n, m;
    vector<vector<T>> tree;

    FenwickPointUpdate2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<T>(m + 1, 0)) {}

    void update(int i, int j, T v){
        for (int x = i; x <= n; x += x & -x){
            for (int y = j; y <= m; y += y & -y){
                tree[x][y] += v;
            }
        }
    }

    T query(int i, int j){
        if (i <= 0 || j <= 0) return 0;
        T res = 0;
        for (int x = i; x > 0; x -= x & -x){
            for (int y = j; y > 0; y -= y & -y){
                res += tree[x][y];
            }
        }
        return res;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        return query(k, l) - query(i - 1, l) - query(k, j - 1) + query(i - 1, j - 1);
    }
};

/***
 * Range updates, point queries
 * update(i, j, k, l, v): Add v to rectangle from [i, j] to [k, l]
 * query(i, j): Return value at index [i, j]
***/

template <typename T>
struct FenwickRangeUpdate2D{
    int n, m;
    vector<vector<T>> tree;

    FenwickRangeUpdate2D(int n = 0, int m = 0) : n(n), m(m), tree(n + 1, vector<T>(m + 1, 0)) {}

    void update(int i, int j, T v){
        for (int x = i; x <= n; x += x & -x){
            for (int y = j; y <= m; y += y & -y){
                tree[x][y] += v;
            }
        }
    }

    void update(int i, int j, int k, int l, T v){
        if (i > k || j > l) return;
        update(i, j, v);
        update(k + 1, j, -v);
        update(k + 1, l + 1, v);
        update(i, l + 1, -v);
    }

    T query(int i, int j){
        T res = 0;
        for (int x = i; x > 0; x -= x & -x){
            for (int y = j; y > 0; y -= y & -y){
                res += tree[x][y];
            }
        }
        return res;
    }
};

/***
 * Range updates, range queries
 * update(i, j, k, l, v): Add v to rectangle from [i, j] to [k, l]
 * query(i, j, k, l): Return sum of rectangle from [i, j] to [k, l]
***/

template <typename T>
struct FenwickFull2D{
    int n, m;
    vector<vector<T>> tree[4];

    FenwickFull2D(int n = 0, int m = 0) : n(n), m(m) {
        for (int i = 0; i < 4; i++){
            tree[i].assign(n + 1, vector<T>(m + 1, 0));
        }
    }

    void update(int p, int q, T v){
        if (p <= 0 || q <= 0 || p > n || q > m) return;

        T c = p - 1, d = q - 1;
        for (int i = p; i <= n; i += i & -i){
            for (int j = q; j <= m; j += j & -j){
                tree[0][i][j] += v;
                tree[1][i][j] += v * d;
                tree[2][i][j] += v * c;
                tree[3][i][j] += v * c * d;
            }
        }
    }

    void update(int i, int j, int k, int l, T v){
        if (i > k || j > l) return;
        update(i, j, v);
        update(k + 1, j, -v);
        update(k + 1, l + 1, v);
        update(i, l + 1, -v);
    }

    T query(int p, int q){
        if (p <= 0 || q <= 0) return 0;

        T x = 0, y = 0, z = 0;
        for (int i = p; i > 0; i -= i & -i){
            T c = 0, d = 0;
            for (int j = q; j > 0; j -= j & -j){
                c += tree[0][i][j];
                d += tree[1][i][j];
                y += tree[2][i][j];
                z += tree[3][i][j];
            }
            x += c * q - d;
        }
        return x * p - y * q + z;
    }

    T query(int i, int j, int k, int l){
        if (i > k || j > l) return 0;
        return query(k, l) - query(i - 1, l) - query(k, j - 1) + query(i - 1, j - 1);
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate2D<int>(10, 10);
    fen.update(2, 3, 5);
    fen.update(5, 7, 3);
    fen.update(4, 4, 7);

    assert(fen.query(2, 3, 2, 3) == 5);
    assert(fen.query(1, 1, 10, 10) == 15);
    assert(fen.query(4, 4, 6, 8) == 10);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate2D<int>(10, 10);
    fen2.update(2, 3, 5, 7, 4);
    fen2.update(4, 4, 8, 9, 2);

    assert(fen2.query(2, 3) == 4);
    assert(fen2.query(5, 6) == 6);
    assert(fen2.query(4, 4) == 6);
    assert(fen2.query(1, 1) == 0);

    /// Range update, range query
    auto fen3 = FenwickFull2D<long long>(10, 10);
    fen3.update(2, 3, 5, 7, 4);
    fen3.update(4, 4, 8, 9, 2);

    assert(fen3.query(2, 3, 5, 7) == 96);
    assert(fen3.query(4, 4, 8, 9) == 92);
    assert(fen3.query(2, 3, 8, 9) == 140);

    return 0;
}
