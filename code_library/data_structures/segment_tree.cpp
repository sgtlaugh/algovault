/***
 *
 * Segment Tree with Lazy Propagation
 * Fast range updates and range queries
 *
 * Complexity: O(log n) per update/query/descent
 *
 * Default: range add update, range sum query, positions are 1-indexed
 * To customize: modify lines marked with // CHANGE
 *
 * Descent (binary search on the tree), pred must be monotone: true on the empty range, and once false it stays false
 * as the range grows
 *   max_right(l, pred): largest r in [l - 1, n] with pred(query(l, r)) true, r = l - 1 meaning the empty range
 *   min_left(r, pred): smallest l in [1, r + 1] with pred(query(l, r)) true, l = r + 1 meaning the empty range
 *   Example: with non-negative values, max_right(l, [&](long long s){ return s <= k; }) is the last r with sum <= k
 *   Descent starts from t_id, so pass the real identity to whichever constructor is used (the default T() is 0,
 *   wrong for min/max: a min tree built without INT_MAX returns wrong positions whenever pred(0) holds)
 *
 * Examples:
 *   Range min with range set: change merge to min(), apply to lz, compose to new_lz, initial value t_id to INT_MAX
 *   Range max with range set: change merge to max(), apply to lz, compose to new_lz, initial value t_id to INT_MIN
 *   Range gcd with range set: change merge to __gcd(), apply to lz, compose to new_lz, initial value t_id to 0
 *   Arithmetic progression add (x, x + d, x + 2d, ... on [l, r]), range sum: lazy {c0, c1} adds c0 + c1 * i at index i
 *     T = {sum, cnt, idx_sum}, built with the vector constructor from {ar[i - 1], 1, i} so leaves know their index,
 *     t_id = {0, 0, 0}
 *     merge and compose add componentwise, apply returns {sum + c0 * cnt + c1 * idx_sum, cnt, idx_sum}
 *     update(l, r, {x - d * l, d})
 *     c1 * idx_sum reaches d * n^2 / 2 (5e18 at d = 1e9, n = 1e5) and composed lazies add up, so large inputs
 *     overflow long long: use a modulus or __int128
 *
 * A pending update is tracked with a flag rather than a reserved lazy value, so any value (including 0) can be set
 *
***/

#include <bits/stdc++.h>

using namespace std;

template<typename T, typename LazyT = T>
struct SegmentTree {
    int n;
    vector<T> tree;
    vector<LazyT> lazy;
    vector<char> has_lazy;
    T t_id;

    /// t_id must be the identity of merge (0 for sum, INT_MAX for min), every node starts as it
    SegmentTree(int n, T t_id = T()) : n(n), tree(n << 2, t_id), lazy(n << 2), has_lazy(n << 2, 0), t_id(t_id) {}

    SegmentTree(const vector<T>& ar, T t_id = T()) : n(ar.size()), tree(n << 2), lazy(n << 2), has_lazy(n << 2, 0),
                                                     t_id(t_id){
        if (n) build(ar, 1, 1, n);
    }

    T merge(T a, T b){
        return a + b;   // CHANGE: a + b (sum), min(a,b) (min), max(a,b) (max), __gcd(a,b) (gcd)
    }

    T apply(T val, LazyT lz, int len){
        return val + lz * len;  // CHANGE: val + lz * len (add), lz (set)
    }

    LazyT compose(LazyT old_lz, LazyT new_lz){
        return old_lz + new_lz; // CHANGE: old_lz + new_lz (add), new_lz (set/overwrite)
    }

    void build(const vector<T>& ar, int idx, int a, int b){
        if (a == b){
            tree[idx] = ar[a - 1];
            return;
        }

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        build(ar, p, a, c);
        build(ar, q, c + 1, b);
        tree[idx] = merge(tree[p], tree[q]);
    }

    void push_lazy(int idx, LazyT val){
        lazy[idx] = has_lazy[idx] ? compose(lazy[idx], val) : val;
        has_lazy[idx] = 1;
    }

    void propagate(int idx, int a, int b){
        if (has_lazy[idx]){
            tree[idx] = apply(tree[idx], lazy[idx], b - a + 1);
            if (a != b){
                push_lazy(idx << 1, lazy[idx]);
                push_lazy(idx << 1 | 1, lazy[idx]);
            }
            has_lazy[idx] = 0;
        }
    }

    void update(int idx, int a, int b, int l, int r, LazyT val){
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
        tree[idx] = merge(tree[p], tree[q]);
    }

    T query(int idx, int a, int b, int l, int r){
        propagate(idx, a, b);
        if (a == l && b == r) return tree[idx];

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c) return query(p, a, c, l, r);
        else if (l > c) return query(q, c + 1, b, l, r);
        else return merge(query(p, a, c, l, c), query(q, c + 1, b, c + 1, r));
    }

    /// first position in [max(a, l), b] where pred(acc merged with values l..pos) fails, or -1, acc absorbs every passing node
    template<typename Pred>
    int first_fail(int idx, int a, int b, int l, Pred& pred, T& acc){
        if (b < l) return -1;
        propagate(idx, a, b);
        if (l <= a){
            T nxt = merge(acc, tree[idx]);
            if (pred(nxt)){
                acc = nxt;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = first_fail(idx << 1, a, c, l, pred, acc);
        return res != -1 ? res : first_fail(idx << 1 | 1, c + 1, b, l, pred, acc);
    }

    /// mirror of first_fail scanning leftwards from r, merge order kept as values pos..r for non-commutative merges
    template<typename Pred>
    int last_fail(int idx, int a, int b, int r, Pred& pred, T& acc){
        if (a > r) return -1;
        propagate(idx, a, b);
        if (b <= r){
            T nxt = merge(tree[idx], acc);
            if (pred(nxt)){
                acc = nxt;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = last_fail(idx << 1 | 1, c + 1, b, r, pred, acc);
        return res != -1 ? res : last_fail(idx << 1, a, c, r, pred, acc);
    }

    void update(int l, int r, LazyT val){ update(1, 1, n, l, r, val); }
    T query(int l, int r){ return query(1, 1, n, l, r); }

    template<typename Pred>
    int max_right(int l, Pred pred){
        assert(1 <= l && l <= n + 1 && pred(t_id));
        T acc = t_id;
        int pos = first_fail(1, 1, n, l, pred, acc);
        return pos == -1 ? n : pos - 1;
    }

    template<typename Pred>
    int min_left(int r, Pred pred){
        assert(0 <= r && r <= n && pred(t_id));
        T acc = t_id;
        int pos = last_fail(1, 1, n, r, pred, acc);
        return pos == -1 ? 1 : pos + 1;
    }
};

int main(){
    // Range add, range sum
    auto seg1 = SegmentTree<long long>(10);
    seg1.update(1, 5, 3);
    seg1.update(3, 7, 2);
    assert(seg1.query(1, 5) == 21);
    assert(seg1.query(3, 7) == 19);

    vector<long long> ar = {1, 2, 3, 4, 5};
    auto seg2 = SegmentTree<long long>(ar);
    assert(seg2.query(1, 5) == 15);
    seg2.update(2, 4, 10);
    assert(seg2.query(1, 5) == 45);

    auto seg3 = SegmentTree<long long>(ar);
    auto at_most = [](long long k){ return [k](long long s){ return s <= k; }; };
    assert(seg3.max_right(1, at_most(6)) == 3);
    assert(seg3.max_right(1, at_most(5)) == 2);
    assert(seg3.max_right(2, at_most(1)) == 1);
    assert(seg3.max_right(1, at_most(100)) == 5);
    assert(seg3.max_right(6, at_most(0)) == 5);
    assert(seg3.min_left(5, at_most(9)) == 4);
    assert(seg3.min_left(5, at_most(4)) == 6);
    assert(seg3.min_left(3, at_most(100)) == 1);
    assert(seg3.min_left(0, at_most(0)) == 1);

    seg3.update(2, 4, 10);
    assert(seg3.max_right(1, at_most(13)) == 2);
    assert(seg3.max_right(3, at_most(26)) == 3);
    assert(seg3.min_left(5, at_most(19)) == 4);
    assert(seg3.min_left(4, at_most(40)) == 1);

    return 0;
}
