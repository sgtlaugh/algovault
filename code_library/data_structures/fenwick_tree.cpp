/***
 *
 * Fenwick Tree (Binary Indexed Tree)
 * Supports both point update + range query and range update + point/range query variants
 * The point update variant also finds the first prefix reaching a given sum (lower_bound)
 * 0-based positions, ranges [l, r] inclusive, a position outside [0, n - 1] asserts
 * l > r is an empty range: updates on it do nothing and query returns 0
 *
 * Complexity:
 *   - O(n) to construct, all zeros or from a vector
 *   - O(log n) per update/query/lower_bound
 *
 * FenwickPointUpdate<long long> f(a);   // a[i] at position i
 * f.update(p, v); f.query(l, r); f.prefix(p);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Point updates, range queries
 * FenwickPointUpdate(n): n zeros, FenwickPointUpdate(a): a[i] at position i
 * update(p, v): Add v to position p
 * prefix(p): Return sum of [0, p], 0 for p = -1
 * query(l, r): Return sum of range [l, r]
 * lower_bound(s): Smallest p with prefix(p) >= s, n if none, 0 if s <= 0
 *   Binary lifting in O(log n), requires all values to be non-negative
***/

template <typename T>
struct FenwickPointUpdate{
    int n;
    vector<T> tree;

    FenwickPointUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    /// O(n): each node pushes its finished sum into its parent
    FenwickPointUpdate(const vector<T>& a) : n(a.size()), tree(n + 1, 0){
        for (int i = 1; i <= n; i++){
            tree[i] += a[i - 1];
            if (i + (i & -i) <= n) tree[i + (i & -i)] += tree[i];
        }
    }

    void update(int p, T v){
        assert(0 <= p && p < n);
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }

    T prefix(int p){
        assert(-1 <= p && p < n);
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

    int lower_bound(T s){
        int p = 0, step = 1;
        while (step <= n / 2) step *= 2;

        for (; step > 0; step >>= 1){
            if (p + step <= n && tree[p + step] < s){
                p += step;
                s -= tree[p];
            }
        }

        return p;
    }
};

/***
 * Range updates, point queries
 * update(l, r, v): Add v to range [l, r]
 * query(p): Return value at position p
***/

template <typename T>
struct FenwickRangeUpdate{
    int n;
    vector<T> tree;

    FenwickRangeUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    void update(int l, int r, T v){
        if (l > r) return;
        assert(0 <= l && r < n);
        add_suffix(l, v);
        add_suffix(r + 1, -v);
    }

    T query(int p){
        assert(0 <= p && p < n);
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

private:
    void add_suffix(int p, T v){
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }
};

/***
 * Range updates, range queries
 * update(l, r, v): Add v to range [l, r]
 * prefix(p): Return sum of [0, p], 0 for p = -1
 * query(l, r): Return sum of range [l, r]
***/

template <typename T>
struct FenwickFull{
    int n;
    vector<T> tree, aux;

    FenwickFull(int n = 0) : n(n), tree(n + 1, 0), aux(n + 1, 0) {}

    void update(int l, int r, T v){
        if (l > r) return;
        assert(0 <= l && r < n);
        update_tree(tree, l, v);
        update_tree(tree, r + 1, -v);
        update_tree(aux, l, v * l);
        update_tree(aux, r + 1, -v * (r + 1));
    }

    T prefix(int p){
        assert(-1 <= p && p < n);
        return query_tree(tree, p) * (p + 1) - query_tree(aux, p);
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

private:
    void update_tree(vector<T>& t, int p, T v){
        for (p++; p <= n; p += p & -p) t[p] += v;
    }

    T query_tree(vector<T>& t, int p){
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += t[p];
        return res;
    }
};

int main(){
    /// Point update, range query
    auto fen = FenwickPointUpdate<int>(10);
    fen.update(0, 5);
    fen.update(2, 7);
    fen.update(4, 3);

    assert(fen.query(0, 0) == 5);
    assert(fen.query(0, 2) == 12);
    assert(fen.query(2, 4) == 10);
    assert(fen.query(0, 9) == 15);
    assert(fen.query(3, 2) == 0);
    assert(fen.prefix(-1) == 0);
    assert(fen.prefix(3) == 12);
    assert(fen.prefix(9) == 15);

    assert(fen.lower_bound(-3) == 0);
    assert(fen.lower_bound(0) == 0);
    assert(fen.lower_bound(1) == 0);
    assert(fen.lower_bound(5) == 0);
    assert(fen.lower_bound(6) == 2);
    assert(fen.lower_bound(12) == 2);
    assert(fen.lower_bound(13) == 4);
    assert(fen.lower_bound(15) == 4);
    assert(fen.lower_bound(16) == 10);

    vector<int> values = {5, 0, 7, 0, 3, 0, 0, 0, 0, 0};
    auto built = FenwickPointUpdate<int>(values);
    for (int l = 0; l < 10; l++){
        for (int r = l; r < 10; r++) assert(built.query(l, r) == fen.query(l, r));
    }
    assert(FenwickPointUpdate<int>(vector<int>{}).n == 0);

    auto zeros = FenwickPointUpdate<long long>(4);
    assert(zeros.lower_bound(0) == 0);
    assert(zeros.lower_bound(1) == 4);

    auto empty = FenwickPointUpdate<int>(0);
    assert(empty.lower_bound(1) == 0);

    /// Range update, point query
    auto fen2 = FenwickRangeUpdate<int>(10);
    fen2.update(0, 4, 3);
    fen2.update(2, 6, 2);
    fen2.update(9, 9, 4);

    assert(fen2.query(0) == 3);
    assert(fen2.query(2) == 5);
    assert(fen2.query(4) == 5);
    assert(fen2.query(6) == 2);
    assert(fen2.query(7) == 0);
    assert(fen2.query(9) == 4);

    /// Range update, range query
    auto fen3 = FenwickFull<long long>(10);
    fen3.update(0, 4, 3);
    fen3.update(2, 6, 2);
    fen3.update(9, 9, 4);

    assert(fen3.query(0, 0) == 3);
    assert(fen3.query(0, 2) == 11);
    assert(fen3.query(2, 4) == 15);
    assert(fen3.query(4, 6) == 9);
    assert(fen3.query(0, 9) == 29);
    assert(fen3.prefix(-1) == 0);
    assert(fen3.prefix(1) == 6);

    return 0;
}
