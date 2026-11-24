/***
 *
 * Historic Segment Tree
 * Range add with range sum and max, plus the largest value each cell ever held and its sum over time
 *
 * Complexity: O(log n) per call, O(n) build and memory
 *
 * HistoricSegmentTree st(n): n zeros, n >= 1
 * HistoricSegmentTree st(a): the values of a, a.size() >= 1
 * st.add(l, r, x): a[i] += x for l <= i <= r
 * st.tick(): h[i] += a[i] for every i, one time step of the historic sum
 * st.query_historic_max(l, r): maximum over l <= i <= r of the largest value a[i] ever held, its first value included
 * st.query_historic_sum(l, r): sum of h[l..r], where h starts at 0 and only tick changes it
 * st.query_max(l, r): maximum of a[l..r]
 * st.query_sum(l, r): sum of a[l..r]
 *
 * Notes:
 *   - Positions are 0-indexed and every call needs 0 <= l <= r < n
 *   - Every value ever held and every range sum, historic ones included, must fit in long long
 *   - Historic minimum: store -a and negate query_historic_max
 *   - Sweep recipe: cell l holds f(a[l..r]); for each r, add r's change, tick, answer the queries ending at r
 *   - query_historic_sum(l, r) is then the sum of f over all subarrays inside [l, r], query_historic_max the best one
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct HistoricSegmentTree{
    /// Pending adds and ticks: total add, highest prefix add, tick count, and the add so far summed over ticks
    struct Tag{
        long long add = 0, peak = 0, ticks = 0, ticked_add = 0;
    };

    struct Node{
        long long sum = 0, max = LLONG_MIN, historic_max = LLONG_MIN, historic_sum = 0;
        int len = 0;
        Tag tag;
    };

    int n;
    vector<Node> tree;

    HistoricSegmentTree(int n) : HistoricSegmentTree(vector<long long>(n, 0)) {}

    HistoricSegmentTree(const vector<long long>& a) : n(a.size()){
        assert(n >= 1);
        tree.assign(2 << __lg(2 * n - 1), Node());
        build(a, 1, 0, n - 1);
    }

    void add(int l, int r, long long x){
        assert(0 <= l && l <= r && r < n);
        update(1, 0, n - 1, l, r, Tag{x, max(x, 0LL), 0, 0});
    }

    void tick(){
        apply(1, Tag{0, 0, 1, 0});
    }

    long long query_historic_max(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r).historic_max;
    }

    long long query_historic_sum(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r).historic_sum;
    }

    long long query_max(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r).max;
    }

    long long query_sum(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r).sum;
    }

    void build(const vector<long long>& a, int idx, int l, int r){
        tree[idx].len = r - l + 1;
        if (l == r){
            tree[idx].sum = tree[idx].max = tree[idx].historic_max = a[l];
            return;
        }

        int mid = (l + r) >> 1;
        build(a, idx << 1, l, mid);
        build(a, idx << 1 | 1, mid + 1, r);
        tree[idx] = merge(tree[idx << 1], tree[idx << 1 | 1], tree[idx].len);
    }

    /// The node's own counters use its values from before t, and t runs after the node's pending tag
    void apply(int idx, const Tag& t){
        Node& v = tree[idx];
        v.historic_sum += t.ticks * v.sum + t.ticked_add * v.len;
        v.historic_max = max(v.historic_max, v.max + t.peak);
        v.sum += t.add * v.len;
        v.max += t.add;

        Tag& p = v.tag;
        p.ticked_add += t.ticked_add + t.ticks * p.add;
        p.peak = max(p.peak, p.add + t.peak);
        p.ticks += t.ticks;
        p.add += t.add;
    }

    void push(int idx){
        apply(idx << 1, tree[idx].tag);
        apply(idx << 1 | 1, tree[idx].tag);
        tree[idx].tag = Tag();
    }

    static Node merge(const Node& a, const Node& b, int len){
        Node res;
        res.sum = a.sum + b.sum;
        res.max = max(a.max, b.max);
        res.historic_max = max(a.historic_max, b.historic_max);
        res.historic_sum = a.historic_sum + b.historic_sum;
        res.len = len;
        return res;
    }

    void update(int idx, int l, int r, int ql, int qr, const Tag& t){
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr){
            apply(idx, t);
            return;
        }

        push(idx);
        int mid = (l + r) >> 1;
        update(idx << 1, l, mid, ql, qr, t);
        update(idx << 1 | 1, mid + 1, r, ql, qr, t);
        tree[idx] = merge(tree[idx << 1], tree[idx << 1 | 1], tree[idx].len);
    }

    Node query(int idx, int l, int r, int ql, int qr){
        if (ql <= l && r <= qr) return tree[idx];

        push(idx);
        int mid = (l + r) >> 1;
        if (qr <= mid) return query(idx << 1, l, mid, ql, qr);
        if (mid < ql) return query(idx << 1 | 1, mid + 1, r, ql, qr);
        return merge(query(idx << 1, l, mid, ql, qr), query(idx << 1 | 1, mid + 1, r, ql, qr), 0);
    }
};

int main(){
    /// a = {3, -1, 4}: add 5 to a[1..2], tick, add -10 to everything, tick
    HistoricSegmentTree st(vector<long long>{3, -1, 4});
    st.add(1, 2, 5);
    assert(st.query_max(0, 2) == 9);

    st.tick();
    st.add(0, 2, -10);
    st.tick();
    assert(st.query_sum(0, 2) == -14);            /// a = {-7, -6, -1}
    assert(st.query_max(0, 2) == -1);
    assert(st.query_historic_max(0, 0) == 3);     /// a[0] went 3, 3, -7
    assert(st.query_historic_max(1, 2) == 9);     /// a[2] went 4, 9, -1
    assert(st.query_historic_sum(0, 2) == 2);     /// first tick {3, 4, 9}, second {-7, -6, -1}: 16 - 14
    assert(st.query_historic_sum(2, 2) == 8);     /// 9 + -1

    return 0;
}
