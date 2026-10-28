/***
 *
 * Fenwick Tree (Binary Indexed Tree)
 * Sums over an array under updates: point add and range sum, range add and point value, or range add and range sum
 *
 * Complexity:
 *   - build: O(n)
 *   - update, prefix, query, lower_bound: O(log n)
 *
 * FenwickPointUpdate<T> fen(n): point add, range sum, over n zeros
 * FenwickPointUpdate<T> fen(a): the same over a copy of a
 * fen.update(p, v): add v to a[p]
 * fen.prefix(p): sum of a[0..p], 0 for p = -1
 * fen.query(l, r): sum of a[l..r], 0 when l > r
 * fen.lower_bound(s): smallest p with prefix(p) >= s, n if none, every a[i] must be non-negative
 *
 * FenwickRangeUpdate<T> fen(n): range add, point value, over n zeros
 * fen.update(l, r, v): add v to a[l..r], nothing when l > r
 * fen.query(p): a[p]
 *
 * FenwickFull<T> fen(n): range add, range sum, over n zeros
 * fen.update(l, r, v): add v to a[l..r], nothing when l > r
 * fen.prefix(p): sum of a[0..p], 0 for p = -1
 * fen.query(l, r): sum of a[l..r], 0 when l > r
 *
 * Positions are 0-based, a position outside [0, n - 1] asserts
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

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

    /// Binary lifting: p stays a multiple of 2 * step, so tree[p + step] is the sum of a[p .. p + step - 1]
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

    void add_suffix(int p, T v){
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }
};

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

    /// With d the difference array, sum of a[0..p] = (p + 1) * sum d[i] - sum d[i] * i; tree holds d[i], aux holds d[i] * i
    T prefix(int p){
        assert(-1 <= p && p < n);
        return query_tree(tree, p) * (p + 1) - query_tree(aux, p);
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

    void update_tree(vector<T>& t, int p, T v){
        for (p++; p <= n; p += p & -p) t[p] += v;
    }

    T query_tree(vector<T>& t, int p){
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += t[p];
        return res;
    }
};

const int N = 1000000;

int main(){
    vector<bool> composite(N + 1, false);
    vector<int> primes;
    for (int i = 2; i <= N; i++){
        if (composite[i]) continue;
        primes.push_back(i);
        for (long long j = (long long)i * i; j <= N; j += i) composite[j] = true;
    }

    vector<int> rev;
    for (int p : primes){
        int d = 0, digits = 0;
        for (; p; p /= 10, digits++) d = d * 10 + p % 10;
        for (; digits < 7; digits++) d *= 10;
        rev.push_back(d);
    }
    sort(rev.begin(), rev.end());

    int n = rev.size();
    vector<int> factors(n);
    for (int i = 0; i < n; i++){
        int x = rev[i], d = 0;
        for (int j = 0; primes[j] * primes[j] <= x; ){
            if (x % primes[j] == 0) x /= primes[j], d++;
            else j++;
        }
        factors[i] = d + (x != 1);
    }

    FenwickPointUpdate<int> alive(vector<int>(n, 1)), sum(factors);
    char type;
    int x;
    while (scanf(" %c %d", &type, &x) == 2){
        if (type == 'd'){
            int p = lower_bound(rev.begin(), rev.end(), x) - rev.begin();
            alive.update(p, -1);
            sum.update(p, -factors[p]);
        }
        else printf("%d\n", sum.prefix(alive.lower_bound(x + 1)));
    }

    return 0;
}
