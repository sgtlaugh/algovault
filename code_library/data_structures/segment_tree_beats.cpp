/***
 *
 * Segment Tree Beats
 * Range chmin, range chmax, range add and range assign, with range sum, min and max queries
 *
 * Complexity: O(log^2 n) amortized per update, O(log n) per query, O(n) build and memory
 *
 * SegmentTreeBeats st(n): n zeros, SegmentTreeBeats st(a): the values of a, n >= 1
 * Positions are 1-indexed and every range [l, r] is inclusive
 *
 * chmin(l, r, x): a[i] = min(a[i], x) for i in [l, r]
 * chmax(l, r, x): a[i] = max(a[i], x) for i in [l, r]
 * add(l, r, x): a[i] += x, assign(l, r, x): a[i] = x
 * query_sum(l, r), query_min(l, r), query_max(l, r)
 *
 * Each node keeps its maximum, its count and the strict second maximum (same for the minimum)
 * A chmin with second max < x < max only rewrites the maximums, so it is tagged without recursing
 *
 * Values: max |a[i]| plus the sum of |x| over all updates must be at most 9e18, and a queried sum must fit in long long
 * Node sums are kept unsigned so they wrap instead of overflowing when an intermediate sum does not fit
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SegmentTreeBeats{
    enum Op{ADD, ASSIGN, CHMAX, CHMIN};

    /// Marks a missing second maximum / second minimum, values never reach these under the documented range
    static constexpr long long NO_MAX = LLONG_MIN, NO_MIN = LLONG_MAX;

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

    int n;
    vector<Node> tree;

    SegmentTreeBeats(int n) : SegmentTreeBeats(vector<long long>(n, 0)) {}

    /// The mid-split recursion on [1, n] never reaches index 2 * next_pow2(n), 4n would waste up to half
    SegmentTreeBeats(const vector<long long>& a) : n(a.size()) {
        assert(n >= 1);
        tree.assign(2 << __lg(2 * n - 1), Node());
        build(a, 1, 1, n);
    }

    void add(int l, int r, long long x){ update(1, 1, n, l, r, ADD, x); }
    void assign(int l, int r, long long x){ update(1, 1, n, l, r, ASSIGN, x); }
    void chmax(int l, int r, long long x){ update(1, 1, n, l, r, CHMAX, x); }
    void chmin(int l, int r, long long x){ update(1, 1, n, l, r, CHMIN, x); }
    long long query_max(int l, int r){ return query(l, r).max; }
    long long query_min(int l, int r){ return query(l, r).min; }
    long long query_sum(int l, int r){ return (long long)query(l, r).sum; }

    void apply(int idx, Op op, long long x){
        if (op == ADD) apply_add(idx, x);
        else if (op == ASSIGN) apply_assign(idx, x);
        else if (op == CHMAX) apply_chmax(idx, x);
        else apply_chmin(idx, x);
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

    /// Requires min1 < x < min2: only the minimums change, and a node with one or two distinct values also moves its max side
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

    void build(const vector<long long>& a, int idx, int l, int r){
        tree[idx].len = r - l + 1;
        if (l == r){
            Node& v = tree[idx];
            v.sum = a[l - 1], v.max1 = v.min1 = a[l - 1], v.cnt_max = v.cnt_min = 1;
            return;
        }

        int c = (l + r) >> 1;
        build(a, idx << 1, l, c);
        build(a, idx << 1 | 1, c + 1, r);
        pull(idx);
    }

    static bool can_tag(const Node& v, Op op, long long x){
        if (op == CHMIN) return v.max2 < x;
        if (op == CHMAX) return v.min2 > x;
        return true;
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

    static bool no_effect(const Node& v, Op op, long long x){
        return (op == CHMIN && v.max1 <= x) || (op == CHMAX && v.min1 >= x);
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

    /// The pending chmin / chmax of a node is implicit: a child whose max exceeds the node's max gets clamped down to it
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

    Summary query(int l, int r){
        assert(1 <= l && l <= r && r <= n);

        Summary res;
        collect(1, 1, n, l, r, res);
        return res;
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
};

int main(){
    SegmentTreeBeats gorgeous(vector<long long>{1, 2, 3, 4, 5});
    assert(gorgeous.query_max(1, 5) == 5);
    assert(gorgeous.query_sum(1, 5) == 15);
    gorgeous.chmin(3, 5, 3);
    assert(gorgeous.query_max(1, 5) == 3);
    assert(gorgeous.query_sum(1, 5) == 12);

    SegmentTreeBeats st(vector<long long>{5, -2, 7, 7, 0, 3});
    st.chmin(1, 6, 4);
    assert(st.query_sum(1, 6) == 13);
    assert(st.query_max(1, 6) == 4);
    st.chmax(2, 5, 1);
    assert(st.query_sum(1, 6) == 17);
    assert(st.query_min(1, 6) == 1);
    st.add(1, 3, -10);
    assert(st.query_sum(1, 6) == -13);
    assert(st.query_min(1, 3) == -9);
    assert(st.query_max(1, 3) == -6);
    st.assign(3, 4, 0);
    assert(st.query_sum(1, 6) == -11);
    st.chmax(1, 6, -7);
    assert(st.query_sum(1, 6) == -9);
    assert(st.query_min(1, 6) == -7);
    st.add(3, 6, 2);
    st.chmin(2, 5, 2);
    assert(st.query_sum(2, 5) == -1);
    assert(st.query_max(1, 6) == 5);

    SegmentTreeBeats zeros(4);
    assert(zeros.query_sum(1, 4) == 0);
    zeros.chmax(2, 3, 0);
    zeros.chmin(2, 3, 0);
    assert(zeros.query_sum(1, 4) == 0);
    zeros.assign(1, 4, 0);
    zeros.add(4, 4, 6);
    assert(zeros.query_max(1, 4) == 6 && zeros.query_min(1, 4) == 0);

    SegmentTreeBeats single(vector<long long>{42});
    single.chmin(1, 1, 50);
    single.chmax(1, 1, 30);
    assert(single.query_sum(1, 1) == 42);
    single.chmin(1, 1, -1);
    assert(single.query_min(1, 1) == -1);

    const long long E18 = 1000000000000000000LL;
    SegmentTreeBeats big(vector<long long>{4 * E18, 4 * E18, -4 * E18});
    assert(big.query_sum(1, 3) == 4 * E18);
    big.add(1, 3, E18);
    assert(big.query_max(1, 3) == 5 * E18);
    assert(big.query_sum(2, 3) == 2 * E18);
    assert(big.query_sum(1, 3) == 7 * E18);
    big.chmin(1, 2, 0);
    assert(big.query_sum(1, 3) == -3 * E18);
    assert(big.query_min(1, 3) == -3 * E18);

    return 0;
}
