/***
 *
 * Segment Tree with Lazy Propagation
 * Range updates and range queries under any merge and lazy update passed as an Ops struct
 *
 * Complexity: O(n) to build, O(log n) per update/query/descent
 *
 * SegmentTree<T, Ops = AddSum<T>, L = T>, 0-based positions, ranges [l, r] inclusive, asserts 0 <= l <= r < n
 * Default: range add update, range sum query
 *     SegmentTree<long long> st(n); st.update(l, r, x); st.query(l, r);
 * SegmentTree(n, v): n copies of v, and v also becomes t_id, SegmentTree(ar, t_id): built from ar
 * t_id is the identity of merge (0 for sum, INT_MAX for min), read only by hld.cpp's query_path, so with hld
 * build a filled tree through the vector constructor and pass the identity
 *
 * Ops holds three static functions, L is the lazy type when it differs from T:
 *   merge(a, b): fold of two adjacent ranges, a on the left
 *   apply(val, lz, len): fold of a range of len elements after the update lz
 *   compose(old_lz, new_lz): one lazy meaning old_lz followed by new_lz
 * Range set, range min:
 *     struct MinSet{
 *         static int merge(int a, int b){ return min(a, b); }
 *         static int apply(int, int lz, int){ return lz; }
 *         static int compose(int, int lz){ return lz; }
 *     };
 *     SegmentTree<int, MinSet> st(a, INT_MAX);
 * Range add, range min, inheriting compose from AddSum:
 *     struct AddMin : AddSum<long long>{
 *         static long long merge(long long a, long long b){ return min(a, b); }
 *         static long long apply(long long val, long long lz, int){ return val + lz; }
 *     };
 * Range max or gcd with range set: MinSet with max or __gcd as merge, t_id INT_MIN or 0
 * Arithmetic progression add (x, x + d, x + 2d, ... on [l, r]), range sum: lazy {c0, c1} adds c0 + c1 * i at index i
 *     T = {sum, cnt, idx_sum}, L = {c0, c1}, built with the vector constructor from {ar[i], 1, i} so leaves know
 *     their index, merge and compose add componentwise, apply returns {sum + c0 * cnt + c1 * idx_sum, cnt, idx_sum}
 *     update(l, r, {x - d * l, d})
 *     c1 * idx_sum reaches d * n^2 / 2 (5e18 at d = 1e9, n = 1e5) and composed lazies add up, so large inputs
 *     overflow long long: use a modulus or __int128
 *
 * Descent (binary search on the tree), pred must be monotone: once false it stays false as the range grows
 *   max_right(l, pred), l in [0, n]: largest r in [l - 1, n - 1] with pred(query(l, r)) true, r = l - 1 when none
 *   min_left(r, pred), r in [-1, n - 1]: smallest l in [0, r + 1] with pred(query(l, r)) true, l = r + 1 when none
 *   pred is never called on the empty range, so the descent does not depend on t_id
 *   Example: with non-negative values, max_right(l, [&](long long s){ return s <= k; }) is the last r with sum <= k
 *
 * A pending update is tracked with a flag rather than a reserved lazy value, so any value (including 0) can be set
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct AddSum{
    static T merge(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
};

/// Inherits Ops so seg.merge(a, b) stays callable from outside, hld.cpp's query_path relies on it
template<typename T, typename Ops = AddSum<T>, typename L = T>
struct SegmentTree : Ops{
    int n;
    T t_id;
    vector<T> tree;
    vector<L> lazy;
    vector<char> has_lazy;

    SegmentTree(int n, T t_id = T()) : SegmentTree(vector<T>(n, t_id), t_id) {}

    SegmentTree(const vector<T>& ar, T t_id = T()) : n(ar.size()), t_id(t_id), tree(4 * n), lazy(4 * n),
                                                     has_lazy(4 * n){
        if (n) build(ar, 1, 0, n - 1);
    }

    void update(int l, int r, L val){
        assert(0 <= l && l <= r && r < n);
        update(1, 0, n - 1, l, r, val);
    }

    T query(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r);
    }

    template<typename Pred>
    int max_right(int l, Pred pred){
        assert(0 <= l && l <= n);
        T acc = t_id;
        bool empty = true;
        int pos = l < n ? first_fail(1, 0, n - 1, l, pred, acc, empty) : -1;
        return pos == -1 ? n - 1 : pos - 1;
    }

    template<typename Pred>
    int min_left(int r, Pred pred){
        assert(-1 <= r && r < n);
        T acc = t_id;
        bool empty = true;
        int pos = r >= 0 ? last_fail(1, 0, n - 1, r, pred, acc, empty) : -1;
        return pos == -1 ? 0 : pos + 1;
    }

private:
    void build(const vector<T>& ar, int idx, int a, int b){
        if (a == b){
            tree[idx] = ar[a];
            return;
        }

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        build(ar, p, a, c);
        build(ar, q, c + 1, b);
        tree[idx] = Ops::merge(tree[p], tree[q]);
    }

    void push_lazy(int idx, const L& val){
        lazy[idx] = has_lazy[idx] ? Ops::compose(lazy[idx], val) : val;
        has_lazy[idx] = 1;
    }

    void propagate(int idx, int a, int b){
        if (!has_lazy[idx]) return;
        tree[idx] = Ops::apply(tree[idx], lazy[idx], b - a + 1);
        if (a != b){
            push_lazy(idx << 1, lazy[idx]);
            push_lazy(idx << 1 | 1, lazy[idx]);
        }
        has_lazy[idx] = 0;
    }

    void update(int idx, int a, int b, int l, int r, const L& val){
        if (a == l && b == r){
            push_lazy(idx, val);
            propagate(idx, a, b);
            return;
        }
        propagate(idx, a, b);

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c){
            propagate(q, c + 1, b);
            update(p, a, c, l, r, val);
        }
        else if (l > c){
            propagate(p, a, c);
            update(q, c + 1, b, l, r, val);
        }
        else{
            update(p, a, c, l, c, val);
            update(q, c + 1, b, c + 1, r, val);
        }
        tree[idx] = Ops::merge(tree[p], tree[q]);
    }

    T query(int idx, int a, int b, int l, int r){
        propagate(idx, a, b);
        if (a == l && b == r) return tree[idx];

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c) return query(p, a, c, l, r);
        if (l > c) return query(q, c + 1, b, l, r);
        return Ops::merge(query(p, a, c, l, c), query(q, c + 1, b, c + 1, r));
    }

    /// first position in [max(a, l), b] where pred(acc merged with values l..pos) fails, or -1, acc absorbs every passing node
    /// acc starts empty rather than at t_id, so a missing or wrong identity cannot move the answer
    template<typename Pred>
    int first_fail(int idx, int a, int b, int l, Pred& pred, T& acc, bool& empty){
        if (b < l) return -1;
        propagate(idx, a, b);
        if (l <= a){
            T nxt = empty ? tree[idx] : Ops::merge(acc, tree[idx]);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = first_fail(idx << 1, a, c, l, pred, acc, empty);
        return res != -1 ? res : first_fail(idx << 1 | 1, c + 1, b, l, pred, acc, empty);
    }

    /// mirror of first_fail scanning leftwards from r, merge order kept as values pos..r for non-commutative merges
    template<typename Pred>
    int last_fail(int idx, int a, int b, int r, Pred& pred, T& acc, bool& empty){
        if (a > r) return -1;
        propagate(idx, a, b);
        if (b <= r){
            T nxt = empty ? tree[idx] : Ops::merge(tree[idx], acc);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = last_fail(idx << 1 | 1, c + 1, b, r, pred, acc, empty);
        return res != -1 ? res : last_fail(idx << 1, a, c, r, pred, acc, empty);
    }
};

/***
 *
 * Coordinate Compression
 * Compress values in input vector in-place
 * Relative order of values are preserved if make_sorted=true, otherwise ids follow first appearance
 *
 * 0 based indexing for compressed values
 * Returns the distinct values in id order, vals[v[i]] is the original value of v[i]
 *
 * Complexity: O(N log N) in both modes, deterministic (no hashing)
 *
 * Example:
 *   auto vals = compress(v);
 *
***/



template <class T>
vector<T> compress(vector<T>& v, bool make_sorted=true){
    int n = v.size();
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());

    vector<T> vals;
    vector<int> ids(n);
    for (int i = 0; i < n; i++){
        if (!i || order[i - 1].first < order[i].first) vals.push_back(order[i].first);
        ids[order[i].second] = vals.size() - 1;
    }

    /// Relabel by scanning in input order instead of hashing, std::hash is the identity on ints and hackable
    if (!make_sorted){
        vector<int> relabel(vals.size(), -1);
        vector<T> first_seen;
        for (int& id: ids){
            if (relabel[id] == -1){
                relabel[id] = first_seen.size();
                first_seen.push_back(vals[id]);
            }
            id = relabel[id];
        }
        vals.swap(first_seen);
    }

    for (int i = 0; i < n; i++) v[i] = ids[i];
    return vals;
}

/// Leaf i stands for the integers [xs[i], xs[i + 1]): every query boundary is a coordinate, so each block is all in or all out
struct Cell{
    int ones, len;
};

/// Lazy: 0 clear, 1 fill, 2 flip, 3 nothing (two flips cancel)
struct FillFlip{
    static Cell merge(Cell a, Cell b){ return {a.ones + b.ones, a.len + b.len}; }
    static Cell apply(Cell c, int lz, int){
        if (lz == 0) return {0, c.len};
        if (lz == 1) return {c.len, c.len};
        if (lz == 2) return {c.len - c.ones, c.len};
        return c;
    }
    static int compose(int old_lz, int new_lz){
        if (new_lz == 3) return old_lz;
        if (new_lz != 2) return new_lz;
        return old_lz == 0 ? 1 : old_lz == 1 ? 0 : old_lz == 2 ? 3 : 2;
    }
};

int main(){
    int n;
    scanf("%d", &n);
    vector<int> t(n);
    vector<long long> l(n), r(n), pts = {1};
    for (int i = 0; i < n; i++){
        scanf("%d %lld %lld", &t[i], &l[i], &r[i]);
        pts.push_back(l[i]), pts.push_back(r[i] + 1);
    }

    vector<long long> ids = pts;
    vector<long long> xs = compress(ids);
    auto id = [&](long long x){ return int(lower_bound(xs.begin(), xs.end(), x) - xs.begin()); };

    int m = xs.size();
    SegmentTree<Cell, FillFlip, int> st(vector<Cell>(m, {0, 1}), {0, 0});
    for (int i = 0; i < n; i++){
        int a = id(l[i]), b = id(r[i] + 1) - 1;
        st.update(a, b, t[i] == 1 ? 1 : t[i] == 2 ? 0 : 2);
        int full = st.max_right(0, [](Cell c){ return c.ones == c.len; });  /// last leaf of the all-in prefix
        printf("%lld\n", xs[full + 1]);
    }
    return 0;
}
