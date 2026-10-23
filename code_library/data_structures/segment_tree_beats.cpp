/***
 *
 * Segment Tree Beats
 * Range chmin, range chmax, range add and range assign, with range sum, min and max queries
 *
 * Complexity: O(log^2 n) amortized per update, O(log n) per query, O(n) build and memory
 *
 * SegmentTreeBeats st(n): n zeros, n >= 1
 * SegmentTreeBeats st(a): the values of a, a.size() >= 1
 * st.add(l, r, x): a[i] += x for l <= i <= r
 * st.assign(l, r, x): a[i] = x for l <= i <= r
 * st.chmax(l, r, x): a[i] = max(a[i], x) for l <= i <= r
 * st.chmin(l, r, x): a[i] = min(a[i], x) for l <= i <= r
 * st.query_max(l, r): maximum of a[l..r]
 * st.query_min(l, r): minimum of a[l..r]
 * st.query_sum(l, r): sum of a[l..r]
 *
 * Notes:
 *   - Positions are 0-indexed and every call needs 0 <= l <= r < n
 *   - max |a[i]| plus the sum of |x| over all updates must be at most 9e18, and a queried sum must fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SegmentTreeBeats{
    /// Marks a missing second maximum / second minimum, values never reach these under the documented range
    static constexpr long long NO_MAX = LLONG_MIN, NO_MIN = LLONG_MAX;

    /// max2 and min2 are strict, never equal to max1 and min1, and sum is unsigned so an intermediate overflow wraps
    struct Node{
        unsigned long long sum = 0;
        long long max1 = 0, max2 = NO_MAX, min1 = 0, min2 = NO_MIN;
        int cnt_max = 0, cnt_min = 0, len = 0;
        long long add = 0, assign = 0;
        bool has_assign = false;
    };

    struct Summary{
        unsigned long long sum = 0;
        long long min = LLONG_MAX, max = LLONG_MIN;
    };

    enum Op{ADD, ASSIGN, CHMAX, CHMIN};

    int n;
    vector<Node> tree;

    SegmentTreeBeats(int n) : SegmentTreeBeats(vector<long long>(n, 0)) {}

    /// The mid-split recursion on [0, n - 1] never reaches index 2 * next_pow2(n)
    SegmentTreeBeats(const vector<long long>& a) : n(a.size()){
        assert(n >= 1);
        tree.assign(2 << __lg(2 * n - 1), Node());
        build(a, 1, 0, n - 1);
    }

    void add(int l, int r, long long x){ update(l, r, ADD, x); }
    void assign(int l, int r, long long x){ update(l, r, ASSIGN, x); }
    void chmax(int l, int r, long long x){ update(l, r, CHMAX, x); }
    void chmin(int l, int r, long long x){ update(l, r, CHMIN, x); }
    long long query_max(int l, int r){ return query(l, r).max; }
    long long query_min(int l, int r){ return query(l, r).min; }
    long long query_sum(int l, int r){ return (long long)query(l, r).sum; }

    void build(const vector<long long>& a, int idx, int l, int r){
        tree[idx].len = r - l + 1;
        if (l == r){
            Node& v = tree[idx];
            v.sum = a[l], v.max1 = v.min1 = a[l], v.cnt_max = v.cnt_min = 1;
            return;
        }

        int c = (l + r) >> 1;
        build(a, idx << 1, l, c);
        build(a, idx << 1 | 1, c + 1, r);
        pull(idx);
    }

    void update(int l, int r, Op op, long long x){
        assert(0 <= l && l <= r && r < n);
        update(1, 0, n - 1, l, r, op, x);
    }

    void update(int idx, int a, int b, int l, int r, Op op, long long x){
        if (r < a || b < l || no_effect(tree[idx], op, x)) return;
        if (l <= a && b <= r && can_tag(tree[idx], op, x)){
            apply(idx, op, x);
            return;
        }

        push(idx);
        int c = (a + b) >> 1;
        update(idx << 1, a, c, l, r, op, x);
        update(idx << 1 | 1, c + 1, b, l, r, op, x);
        pull(idx);
    }

    Summary query(int l, int r){
        assert(0 <= l && l <= r && r < n);

        Summary res;
        collect(1, 0, n - 1, l, r, res);
        return res;
    }

    void pull(int idx){
        const Node &p = tree[idx << 1], &q = tree[idx << 1 | 1];
        Node& v = tree[idx];
        v.sum = p.sum + q.sum;

        if (p.max1 == q.max1) v.max1 = p.max1, v.cnt_max = p.cnt_max + q.cnt_max, v.max2 = max(p.max2, q.max2);
        else if (p.max1 > q.max1) v.max1 = p.max1, v.cnt_max = p.cnt_max, v.max2 = max(p.max2, q.max1);
        else v.max1 = q.max1, v.cnt_max = q.cnt_max, v.max2 = max(p.max1, q.max2);

        if (p.min1 == q.min1) v.min1 = p.min1, v.cnt_min = p.cnt_min + q.cnt_min, v.min2 = min(p.min2, q.min2);
        else if (p.min1 < q.min1) v.min1 = p.min1, v.cnt_min = p.cnt_min, v.min2 = min(p.min2, q.min1);
        else v.min1 = q.min1, v.cnt_min = q.cnt_min, v.min2 = min(p.min1, q.min2);
    }

    static bool no_effect(const Node& v, Op op, long long x){
        return (op == CHMIN && v.max1 <= x) || (op == CHMAX && v.min1 >= x);
    }

    static bool can_tag(const Node& v, Op op, long long x){
        if (op == CHMIN) return v.max2 < x;
        if (op == CHMAX) return v.min2 > x;
        return true;
    }

    void apply(int idx, Op op, long long x){
        if (op == ADD) apply_add(idx, x);
        else if (op == ASSIGN) apply_assign(idx, x);
        else if (op == CHMAX) apply_chmax(idx, x);
        else apply_chmin(idx, x);
    }

    /// A node's pending chmin / chmax is implicit: a child whose max exceeds the node's max is clamped down to it
    void push(int idx){
        Node& v = tree[idx];
        for (int child : {idx << 1, idx << 1 | 1}){
            if (v.has_assign){
                apply_assign(child, v.assign);
                continue;
            }

            if (v.add) apply_add(child, v.add);
            if (tree[child].max1 > v.max1) apply_chmin(child, v.max1);
            if (tree[child].min1 < v.min1) apply_chmax(child, v.min1);
        }
        v.add = 0, v.has_assign = false;
    }

    void collect(int idx, int a, int b, int l, int r, Summary& res){
        if (r < a || b < l) return;
        if (l <= a && b <= r){
            res.sum += tree[idx].sum;
            res.min = min(res.min, tree[idx].min1), res.max = max(res.max, tree[idx].max1);
            return;
        }

        push(idx);
        int c = (a + b) >> 1;
        collect(idx << 1, a, c, l, r, res);
        collect(idx << 1 | 1, c + 1, b, l, r, res);
    }

    void apply_add(int idx, long long x){
        Node& v = tree[idx];
        v.sum += (unsigned long long)x * v.len;
        v.max1 += x, v.min1 += x;
        if (v.max2 != NO_MAX) v.max2 += x;
        if (v.min2 != NO_MIN) v.min2 += x;

        if (v.has_assign) v.assign += x;
        else v.add += x;
    }

    void apply_assign(int idx, long long x){
        Node& v = tree[idx];
        v.sum = (unsigned long long)x * v.len;
        v.max1 = v.min1 = x, v.max2 = NO_MAX, v.min2 = NO_MIN;
        v.cnt_max = v.cnt_min = v.len;
        v.assign = x, v.has_assign = true, v.add = 0;
    }

    /// Requires min1 < x < min2: only minimums change, a node with one or two distinct values also moves its max side
    void apply_chmax(int idx, long long x){
        Node& v = tree[idx];
        v.sum += ((unsigned long long)x - (unsigned long long)v.min1) * v.cnt_min;
        if (v.max1 == v.min1) v.max1 = x;
        else if (v.max2 == v.min1) v.max2 = x;
        v.min1 = x;

        if (v.has_assign) v.assign = x;
    }

    /// Requires max2 < x < max1, the mirror of apply_chmax
    void apply_chmin(int idx, long long x){
        Node& v = tree[idx];
        v.sum += ((unsigned long long)x - (unsigned long long)v.max1) * v.cnt_max;
        if (v.min1 == v.max1) v.min1 = x;
        else if (v.min2 == v.max1) v.min2 = x;
        v.max1 = x;

        if (v.has_assign) v.assign = x;
    }
};

int main(){
    /// Each update, then queries on the array it leaves
    SegmentTreeBeats st(vector<long long>{5, -2, 7, 7, 0, 3});
    st.chmin(0, 5, 4);                  /// 4 -2 4 4 0 3
    assert(st.query_max(0, 5) == 4);
    assert(st.query_sum(0, 5) == 13);
    st.chmax(1, 4, 1);                  /// 4 1 4 4 1 3
    assert(st.query_min(0, 5) == 1);
    st.add(0, 2, -10);                  /// -6 -9 -6 4 1 3
    assert(st.query_sum(0, 2) == -21);
    st.assign(2, 3, 0);                 /// -6 -9 0 0 1 3
    assert(st.query_sum(0, 5) == -11);
    assert(st.query_max(0, 3) == 0);

    return 0;
}
