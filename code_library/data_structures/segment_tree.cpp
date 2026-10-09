/***
 *
 * Segment Tree with Lazy Propagation
 * Fast range updates and range queries
 *
 * Complexity: O(log n) per update/query
 *
 * Default: range add update, range sum query
 * To customize: modify lines marked with // CHANGE
 *
 * Examples:
 *   Range min with range set: change merge to min(), apply to lz, compose to new_lz, initial value t_id to INT_MAX
 *   Range max with range set: change merge to max(), apply to lz, compose to new_lz, initial value t_id to INT_MIN
 *   Range gcd with range set: change merge to __gcd(), apply to lz, compose to new_lz, initial value t_id to 0
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

    /// t_id must be the identity of merge (0 for sum, INT_MAX for min), every node starts as it
    SegmentTree(int n, T t_id = 0) : n(n), tree(n << 2, t_id), lazy(n << 2), has_lazy(n << 2, 0) {}

    SegmentTree(const vector<T>& ar) : n(ar.size()), tree(n << 2), lazy(n << 2), has_lazy(n << 2, 0) {
        build(ar, 1, 1, n);
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

    void update(int l, int r, LazyT val){ update(1, 1, n, l, r, val); }
    T query(int l, int r){ return query(1, 1, n, l, r); }
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

    return 0;
}
