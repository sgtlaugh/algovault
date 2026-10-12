/***
 *
 * Segment Tree with Lazy Propagation
 * Range update, range query and binary search for any associative merge and lazy update
 *
 * Complexity: O(n) to build, O(log n) per update, query, max_right and min_left
 *
 * SegmentTree<T, Ops = AddSum<T>, L = T>: 0-based, ranges [l, r] inclusive, default range add and range sum
 *   SegmentTree<long long> st(n); st.update(l, r, x); st.query(l, r);
 *   SegmentTree(n, v): n copies of v, SegmentTree(ar, t_id): built from ar, t_id the identity of merge
 *   max_right(l, pred), l in [0, n]: largest r with pred(query(l, r)) true, l - 1 if none
 *   min_left(r, pred), r in [-1, n - 1]: smallest l with pred(query(l, r)) true, r + 1 if none
 *   pred must be monotone: once false it stays false as the range grows
 *
 * Ops holds three static functions, L is the lazy type when it differs from T:
 *   merge(a, b): fold of two adjacent ranges, a on the left
 *   apply(val, lz, len): value of a range of len elements after lz
 *   compose(old, new): one lazy meaning old then new, a lazy that changes nothing must return old
 *
 * Range set, range min:
 *     struct MinSet{
 *         static int merge(int a, int b){ return min(a, b); }
 *         static int apply(int, int lz, int){ return lz; }
 *         static int compose(int, int lz){ return lz; }
 *     };
 *     SegmentTree<int, MinSet> st(a, INT_MAX);
 *
 * Arithmetic progression add (x, x + d, ... on [l, r]) with range sum: T = {sum, cnt, idx_sum}, L = {c0, c1}
 *     leaves {ar[i], 1, i}, merge and compose add componentwise, apply adds c0 * cnt + c1 * idx_sum,
 *     update(l, r, {x - d * l, d}), sums reach d * n^2 / 2 so large inputs need a modulus or __int128
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

int main(){
    SegmentTree<long long> st(vector<long long>{1, 2, 3, 4, 5});
    st.update(1, 3, 10);                                                /// 1 12 13 14 5
    assert(st.query(0, 4) == 45);
    assert(st.query(2, 3) == 27);
    assert(st.max_right(0, [](long long s){ return s <= 13; }) == 1);  /// 1 + 12 fits, adding 13 does not
    assert(st.min_left(4, [](long long s){ return s <= 19; }) == 3);   /// 14 + 5

    struct MinSet{
        static int merge(int a, int b){ return min(a, b); }
        static int apply(int, int lz, int){ return lz; }
        static int compose(int, int lz){ return lz; }
    };
    SegmentTree<int, MinSet> mn(vector<int>{5, 3, 8, 6, 7}, INT_MAX);
    mn.update(0, 2, 9);                                                 /// 9 9 9 6 7
    assert(mn.query(0, 4) == 6);
    assert(mn.query(0, 2) == 9);
    return 0;
}
