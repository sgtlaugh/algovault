/***
 *
 * Li Chao Tree
 * Maintains lines y = k * x + b over integer x in [lo, hi] and answers the minimum (or maximum) at a point
 *
 * Complexity:
 *   - O(log C) per add_line and query, O(log^2 C) per add_segment, where C = hi - lo + 1
 *   - Nodes are created on demand: at most 1 per add_line, O(log C) per add_segment
 *   - add_segment costs ~1 KB per call at C = 2e9, so 2e5 segments over that range take ~200 MB
 *
 * LiChaoTree<T> tree(lo, hi): empty tree over x in [lo, hi], the range can span up to [-1e18, 1e18]
 * add_line(k, b): adds y = k * x + b for every x in [lo, hi]
 * add_segment(k, b, l, r): adds the line only for x in [l, r] inclusive, the part outside [lo, hi] is ignored
 * query(x): minimum over the lines covering x, or NONE if no line covers it
 *
 * LiChaoTree<T, true> answers the maximum instead, and NONE becomes the lowest value of T
 * Every k * x + b over the range must fit in T, e.g. |k|, |x| <= 1e9 and |b| <= 1e18 for long long
 *
 * Convex hull trick: dp[i] = min over j < i of (dp[j] + a[j] * x[i]) becomes
 * dp[i] = tree.query(x[i]), then tree.add_line(a[i], dp[i])
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, bool MAXIMIZE = false>
struct LiChaoTree{
    static constexpr T NONE = MAXIMIZE ? numeric_limits<T>::lowest() : numeric_limits<T>::max();

    struct Line{
        T k, b;

        T eval(long long x) const{
            return k * x + b;
        }
    };

    struct Node{
        Line line;
        bool has;
        int left, right;
    };

    long long lo, hi;
    int root = -1;
    vector<Node> nodes;

    LiChaoTree(long long lo, long long hi) : lo(lo), hi(hi) {}

    void add_line(T k, T b){
        root = insert(root, lo, hi, {k, b});
    }

    void add_segment(T k, T b, long long l, long long r){
        l = max(l, lo), r = min(r, hi);
        if (l <= r) root = insert_segment(root, lo, hi, l, r, {k, b});
    }

    T query(long long x) const{
        assert(lo <= x && x <= hi);

        T res = NONE;
        long long l = lo, r = hi;
        for (int node = root; node != -1;){
            const Node& cur = nodes[node];
            if (cur.has){
                T y = cur.line.eval(x);
                if (better(y, res)) res = y;
            }

            long long m = l + (r - l) / 2;
            if (x <= m) node = cur.left, r = m;
            else node = cur.right, l = m + 1;
        }
        return res;
    }

    static bool better(T a, T b){
        return MAXIMIZE ? a > b : a < b;
    }

    int new_node(const Line& line, bool has){
        nodes.push_back({line, has, -1, -1});
        return (int)nodes.size() - 1;
    }

    /// Walks a single root to leaf path, each node keeps the line that wins at its midpoint
    int insert(int node, long long l, long long r, Line line){
        if (node == -1) return new_node(line, true);

        for (int top = node; ; ){
            if (!nodes[node].has){
                nodes[node].line = line, nodes[node].has = true;
                return top;
            }

            long long m = l + (r - l) / 2;
            Line& cur = nodes[node].line;
            bool left_better = better(line.eval(l), cur.eval(l));
            bool mid_better = better(line.eval(m), cur.eval(m));
            if (mid_better) swap(cur, line);
            if (l == r) return top;

            bool go_left = left_better != mid_better;
            int child = go_left ? nodes[node].left : nodes[node].right;
            if (child == -1){
                child = new_node(line, true);
                (go_left ? nodes[node].left : nodes[node].right) = child;
                return top;
            }

            node = child;
            if (go_left) r = m;
            else l = m + 1;
        }
    }

    int insert_segment(int node, long long l, long long r, long long ql, long long qr, const Line& line){
        if (qr < l || r < ql) return node;
        if (ql <= l && r <= qr) return insert(node, l, r, line);
        if (node == -1) node = new_node(line, false);

        long long m = l + (r - l) / 2;
        int left = insert_segment(nodes[node].left, l, m, ql, qr, line);  /// the call can reallocate nodes, assign after it
        nodes[node].left = left;
        int right = insert_segment(nodes[node].right, m + 1, r, ql, qr, line);
        nodes[node].right = right;
        return node;
    }
};

int main(){
    LiChaoTree<long long> tree(-10, 10);
    assert(tree.query(0) == LiChaoTree<long long>::NONE);

    tree.add_line(1, 0);
    tree.add_line(-1, 0);
    assert(tree.query(3) == -3);
    assert(tree.query(-4) == -4);
    assert(tree.query(0) == 0);

    tree.add_line(0, -2);
    assert(tree.query(1) == -2);
    assert(tree.query(5) == -5);
    assert(tree.query(-10) == -10);

    LiChaoTree<long long> seg(0, 10);
    seg.add_segment(0, 5, 2, 4);
    assert(seg.query(1) == LiChaoTree<long long>::NONE);
    assert(seg.query(2) == 5);
    assert(seg.query(4) == 5);
    assert(seg.query(5) == LiChaoTree<long long>::NONE);

    seg.add_line(1, 0);
    assert(seg.query(3) == 3);
    assert(seg.query(1) == 1);
    assert(seg.query(10) == 10);

    seg.add_segment(-100, 0, 6, 3);
    seg.add_segment(0, -7, -5, 0);
    assert(seg.query(0) == -7);
    assert(seg.query(1) == 1);

    LiChaoTree<long long, true> best(-10, 10);
    assert(best.query(7) == numeric_limits<long long>::lowest());
    best.add_line(2, 1);
    best.add_line(-1, 4);
    assert(best.query(0) == 4);
    assert(best.query(2) == 5);
    assert(best.query(-10) == 14);

    LiChaoTree<long long> point(7, 7);
    point.add_line(3, -1);
    point.add_line(2, 5);
    assert(point.query(7) == 19);

    const long long E18 = 1000000000000000000LL;
    LiChaoTree<long long> wide(-E18, E18);
    wide.add_line(2, 0);
    wide.add_line(-3, E18);
    assert(wide.query(0) == 0);
    assert(wide.query(E18) == -2 * E18);
    assert(wide.query(-E18) == -2 * E18);

    LiChaoTree<double> real(0, 10);
    real.add_line(0.5, 0.25);
    real.add_line(-0.5, 4.0);
    assert(abs(real.query(3) - 1.75) < 1e-9);
    assert(abs(real.query(10) - (-1.0)) < 1e-9);

    return 0;
}
