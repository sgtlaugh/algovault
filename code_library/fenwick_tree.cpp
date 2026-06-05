/***
 *
 * Fenwick Tree (Binary Indexed Tree)
 * Supports both point update + range query and range update + point/range query variants
 * All operations in O(log n)
 * 1-based indexing for elements
 *
 * Complexity:
 *   - O(n) to build
 *   - O(log n) per update/query
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Point updates, range queries
 * update(p, v): Add v to index p
 * query(l, r): Return sum of range [l, r]
***/

template <typename T>
struct FenwickPointUpdate{
    int n;
    vector<T> tree;

    FenwickPointUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    void update(int p, T v){
        for (; p <= n; p += p & -p) tree[p] += v;
    }

    T query(int p){
        T res = 0;
        for (; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

    T query(int l, int r){
        if (l > r) return 0;
        return query(r) - query(l - 1);
    }
};

/***
 * Range updates, point queries
 * update(l, r, v): Add v to range [l, r]
 * query(p): Return value at index p
***/

template <typename T>
struct FenwickRangeUpdate{
    int n;
    vector<T> tree;

    FenwickRangeUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    void update(int p, T v){
        for (; p <= n; p += p & -p) tree[p] += v;
    }

    void update(int l, int r, T v){
        if (l > r) return;
        update(l, v);
        update(r + 1, -v);
    }

    T query(int p){
        T res = 0;
        for (; p > 0; p -= p & -p) res += tree[p];
        return res;
    }
};

/***
 * Range updates, range queries
 * update(l, r, v): Add v to range [l, r]
 * query(l, r): Return sum of range [l, r]
***/

template <typename T>
struct FenwickFull{
    int n;
    vector<T> tree, aux;

    FenwickFull(int n = 0) : n(n), tree(n + 1, 0), aux(n + 1, 0) {}

    void update_tree(vector<T>& t, int p, T v){
        for (; p <= n; p += p & -p) t[p] += v;
    }

    void update(int l, int r, T v){
        if (l > r) return;
        update_tree(tree, l, v);
        update_tree(tree, r + 1, -v);
        update_tree(aux, l, v * (l - 1));
        update_tree(aux, r + 1, -v * r);
    }

    T query_tree(vector<T>& t, int p){
        T res = 0;
        for (; p > 0; p -= p & -p) res += t[p];
        return res;
    }

    T query(int p){
        return query_tree(tree, p) * p - query_tree(aux, p);
    }

    T query(int l, int r){
        if (l > r) return 0;
        return query(r) - query(l - 1);
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate<int>(10);
    fen.update(1, 5);
    fen.update(3, 7);
    fen.update(5, 3);

    assert(fen.query(1, 1) == 5);
    assert(fen.query(1, 3) == 12);
    assert(fen.query(3, 5) == 10);
    assert(fen.query(1, 10) == 15);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate<int>(10);
    fen2.update(1, 5, 3);
    fen2.update(3, 7, 2);

    assert(fen2.query(1) == 3);
    assert(fen2.query(3) == 5);
    assert(fen2.query(5) == 5);
    assert(fen2.query(7) == 2);
    assert(fen2.query(8) == 0);

    /// Range update, range query
    auto fen3 = FenwickFull<long long>(10);
    fen3.update(1, 5, 3);
    fen3.update(3, 7, 2);

    assert(fen3.query(1, 1) == 3);
    assert(fen3.query(1, 3) == 11);
    assert(fen3.query(3, 5) == 15);
    assert(fen3.query(5, 7) == 9);
    assert(fen3.query(1, 10) == 25);

    return 0;
}
