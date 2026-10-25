/***
 *
 * Li Chao Tree
 * Holds lines y = k * x + b over integer x in [lo, hi] and returns the minimum or maximum at a point
 *
 * Complexity:
 *   - LiChaoTree: O(log C) per add_line and query, O(log^2 C) per add_segment, C = hi - lo + 1
 *   - PersistentLiChaoTree: O(log C) per add_line and query, add_line makes at most log2(C) + 2 nodes
 *
 * LiChaoTree<T> tree(lo, hi): no lines yet, x in [lo, hi], the range can span up to [-1e18, 1e18]
 * LiChaoTree<T, true> tree(lo, hi): the same tree answering the maximum
 * tree.add_line(k, b): add y = k * x + b for every x in [lo, hi]
 * tree.add_segment(k, b, l, r): add the line only for x in [l, r] inclusive, the part outside [lo, hi] is ignored
 * tree.query(x): minimum over the lines covering x, or NONE if no line covers it, lo <= x <= hi
 *
 * PersistentLiChaoTree<T> pt(lo, hi): version 0 holds no lines, same range limits as LiChaoTree
 * PersistentLiChaoTree<T, true> pt(lo, hi): the same tree answering the maximum
 * pt.add_line(version, k, b): a new version with the lines of version plus y = k * x + b, older versions stay intact, version is 0 or a value add_line returned
 * pt.query(version, x): minimum over the lines of version at x, or NONE if it has no line, lo <= x <= hi
 *
 * NONE is the largest value of T, or the lowest when answering the maximum
 * Every k * x + b over the range must fit in T, e.g. |k|, |x| <= 1e9 and |b| <= 1e18 for long long
 * LiChaoTree makes nodes on demand, at most 1 per add_line and O(log C) per add_segment
 * A node is 24 bytes for long long, random segments at C = 2e9 average ~31 nodes, so 4e5 of them take ~400 MB peak
 * Convex hull trick: dp[i] = min over j < i of (dp[j] + a[j] * x[i]) is dp[i] = tree.query(x[i]), then tree.add_line(a[i], dp[i])
 * The same on root to node paths: dp[v] = pt.query(ver[parent], x[v]), then ver[v] = pt.add_line(ver[parent], a[v], dp[v])
 * The persistent DFS needs no rollback when it backtracks, siblings branch off the same parent version
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

    /// A node on a segment's split path holds y = NONE, which never wins, so no node needs a has-line flag
    struct Node{
        Line line;
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
            T y = cur.line.eval(x);
            if (better(y, res)) res = y;

            long long m = l + (r - l) / 2;
            if (x <= m) node = cur.left, r = m;
            else node = cur.right, l = m + 1;
        }
        return res;
    }

    static bool better(T a, T b){
        return MAXIMIZE ? a > b : a < b;
    }

    int new_node(const Line& line){
        nodes.push_back({line, -1, -1});
        return (int)nodes.size() - 1;
    }

    /// Walks a single root to leaf path, each node keeps the line that wins at its midpoint
    /// Two lines cross at most once, so when one wins at both ends of a node's range it wins on all of it and the walk stops
    int insert(int node, long long l, long long r, Line line){
        if (node == -1) return new_node(line);

        for (int top = node; ; ){
            Line& cur = nodes[node].line;
            bool left_better = better(line.eval(l), cur.eval(l));
            bool right_better = better(line.eval(r), cur.eval(r));
            if (left_better == right_better){
                if (left_better) cur = line;
                return top;
            }

            long long m = l + (r - l) / 2;
            bool mid_better = better(line.eval(m), cur.eval(m));
            if (mid_better) swap(cur, line);

            bool go_left = left_better != mid_better;
            int child = go_left ? nodes[node].left : nodes[node].right;
            if (child == -1){
                child = new_node(line);
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
        if (node == -1) node = new_node({0, NONE});

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
    using Node = typename Base::Node;

    static constexpr T NONE = Base::NONE;

    long long lo, hi;
    vector<int> roots = {-1};
    vector<Node> nodes;

    PersistentLiChaoTree(long long lo, long long hi) : lo(lo), hi(hi) {}

    /// Midpoint rule of LiChaoTree::insert without its early stop, writing into a copy of every node it visits
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
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cin >> n >> m;
    LiChaoTree<long long, true> tree(0, m);
    for (int i = 0; i < n; i++){
        long long y1, y2;
        cin >> y1 >> y2;
        tree.add_line((y2 - y1) / m, y1);
    }

    for (int x = 0; x <= m; x++) cout << tree.query(x) << (x == m ? '\n' : ' ');
    return 0;
}
