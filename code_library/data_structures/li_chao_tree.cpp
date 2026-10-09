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
 * Persistent Li Chao Tree
 * Every add_line copies its root to leaf path and returns a new version, older versions stay intact
 *
 * Complexity:
 *   - O(log C) per add_line and query
 *   - add_line creates at most log2(C) + 2 nodes
 *
 * PersistentLiChaoTree<T> tree(lo, hi): version 0 is the empty set of lines, same range limits as above
 * add_line(version, k, b): returns a new version holding the lines of version plus y = k * x + b
 * query(version, x): minimum over the lines of version at x, or NONE if it has no line
 * PersistentLiChaoTree<T, true> answers the maximum, same overflow rule as LiChaoTree
 *
 * DP on root to node paths: dp[v] = min over ancestors u of (dp[u] + a[u] * x[v]) becomes
 * dp[v] = tree.query(ver[parent], x[v]), then ver[v] = tree.add_line(ver[parent], a[v], dp[v])
 * No rollback is needed when the DFS backtracks, siblings branch off the same parent version
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

template <typename T, bool MAXIMIZE = false>
struct PersistentLiChaoTree{
    using Base = LiChaoTree<T, MAXIMIZE>;
    using Line = typename Base::Line;

    static constexpr T NONE = Base::NONE;

    struct Node{
        Line line;
        int left, right;
    };

    long long lo, hi;
    vector<int> roots = {-1};
    vector<Node> nodes;

    PersistentLiChaoTree(long long lo, long long hi) : lo(lo), hi(hi) {}

    /// Same descent as LiChaoTree::insert, but writes into a copy of every node it visits
    int add_line(int version, T k, T b){
        assert(0 <= version && version < (int)roots.size());

        Line line = {k, b};
        int node = roots[version], top = -1, parent = -1;
        bool from_left = false;
        long long l = lo, r = hi;

        while (true){
            Node copy = node == -1 ? Node{line, -1, -1} : nodes[node];
            nodes.push_back(copy);
            int cur = (int)nodes.size() - 1;
            if (parent == -1) top = cur;
            else (from_left ? nodes[parent].left : nodes[parent].right) = cur;
            if (node == -1) break;

            long long m = l + (r - l) / 2;
            Line& held = nodes[cur].line;
            bool left_better = Base::better(line.eval(l), held.eval(l));
            bool mid_better = Base::better(line.eval(m), held.eval(m));
            if (mid_better) swap(held, line);
            if (l == r) break;

            parent = cur, from_left = left_better != mid_better;
            node = from_left ? nodes[cur].left : nodes[cur].right;
            if (from_left) r = m;
            else l = m + 1;
        }

        roots.push_back(top);
        return (int)roots.size() - 1;
    }

    T query(int version, long long x) const{
        assert(0 <= version && version < (int)roots.size());
        assert(lo <= x && x <= hi);

        T res = NONE;
        long long l = lo, r = hi;
        for (int node = roots[version]; node != -1;){
            const Node& cur = nodes[node];
            T y = cur.line.eval(x);
            if (Base::better(y, res)) res = y;

            long long m = l + (r - l) / 2;
            if (x <= m) node = cur.left, r = m;
            else node = cur.right, l = m + 1;
        }
        return res;
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

    PersistentLiChaoTree<long long> pt(-10, 10);
    assert(pt.query(0, 4) == PersistentLiChaoTree<long long>::NONE);

    int v1 = pt.add_line(0, 1, 0);
    int v2 = pt.add_line(v1, -1, 0);
    int v3 = pt.add_line(v1, 0, -2);
    int v4 = pt.add_line(v2, 0, -7);
    assert(pt.query(0, 4) == PersistentLiChaoTree<long long>::NONE);
    assert(pt.query(v1, 3) == 3);
    assert(pt.query(v2, 3) == -3);
    assert(pt.query(v2, -4) == -4);
    assert(pt.query(v2, 0) == 0);
    assert(pt.query(v3, 3) == -2);
    assert(pt.query(v3, -5) == -5);
    assert(pt.query(v3, 0) == -2);
    assert(pt.query(v4, 0) == -7);
    assert(pt.query(v4, 10) == -10);
    assert(pt.query(v4, -9) == -9);
    int v5 = pt.add_line(0, 0, 0);
    assert(v5 == 5);
    assert(pt.query(v5, -10) == 0);

    PersistentLiChaoTree<long long, true> pbest(0, 5);
    int p1 = pbest.add_line(0, 2, 1);
    int p2 = pbest.add_line(p1, -1, 6);
    int p3 = pbest.add_line(0, -1, 6);
    assert(pbest.query(0, 2) == numeric_limits<long long>::lowest());
    assert(pbest.query(p1, 0) == 1);
    assert(pbest.query(p2, 0) == 6);
    assert(pbest.query(p2, 1) == 5);
    assert(pbest.query(p2, 5) == 11);
    assert(pbest.query(p3, 5) == 1);

    PersistentLiChaoTree<long long> ppoint(7, 7);
    int q1 = ppoint.add_line(0, 3, -1);
    int q2 = ppoint.add_line(q1, 2, 5);
    assert(ppoint.query(q1, 7) == 20);
    assert(ppoint.query(q2, 7) == 19);

    PersistentLiChaoTree<long long> pwide(-E18, E18);
    int w1 = pwide.add_line(0, 2, 0);
    int w2 = pwide.add_line(w1, -3, E18);
    assert(pwide.query(w1, E18) == 2 * E18);
    assert(pwide.query(w2, E18) == -2 * E18);
    assert(pwide.query(w2, -E18) == -2 * E18);
    assert(pwide.query(w2, 0) == 0);

    return 0;
}
