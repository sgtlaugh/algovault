/***
 *
 * KD-Tree
 * Static 2D point set: the point nearest to a query by Euclidean or Manhattan distance, and the points in a rectangle
 *
 * Complexity:
 *   - KDTree: O(n log n)
 *   - nearest: about O(log n) on random or clustered points, O(n) for a query at the centre of points on a circle
 *   - report: about O(sqrt n + k) for k reported points
 *
 * KDTree tree(points): tree over points whose nearest uses the squared Euclidean distance dx^2 + dy^2
 * KDTree tree(points, KDTree::MANHATTAN): tree over points whose nearest uses the Manhattan distance |dx| + |dy|
 * tree.nearest(q): {dist, i}, a point points[i] nearest to q and its distance, {LLONG_MAX, -1} if there is none
 * tree.nearest(q, skip): {dist, i}, a point nearest to q among all points except points[skip], i = -1 if there is none
 * tree.report(x1, y1, x2, y2): indices of the points with x1 <= x <= x2 and y1 <= y <= y2, in no particular order
 * tree.report(x1, y1, x2, y2, f): calls f(i) for each point with x1 <= x <= x2 and y1 <= y <= y2, building no vector
 * tree.nodes: every point once in tree order, as node.p = points[node.i]
 *
 * Notes:
 *   - |x|, |y| <= 1e9 for points and queries, so a squared distance (at most 8e18) fits in a long long
 *   - nearest(points[i], i) is the nearest other point, at dist 0 when points[i] has a duplicate
 *   - To query at every input point, loop over tree.nodes: tree order is about 2x faster than input order
 *   - Ties between equally near points are broken arbitrarily
 *   - Reimplemented from the idea of the KACTL KD-tree (no licence stated): bounding box pruning, nearer child first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;
};

struct KDTree{
    enum Metric{EUCLIDEAN, MANHATTAN};

    struct Neighbor{
        long long dist;
        int i;
    };

    struct Node{
        Point p;
        int i;
    };

    /// The box of nodes[l, r); left < 0 marks a leaf, else the cells left and right split nodes[l, r) in two
    struct Cell{
        long long x1, x2, y1, y2;
        int l, r, left, right;
    };

    /// Scanning a leaf of 32 points is faster than splitting it further, for nearest and report alike
    static const int LEAF = 32;

    Metric metric;
    vector<Node> nodes;
    vector<Cell> cells;

    KDTree(const vector<Point>& points, Metric metric = EUCLIDEAN) : metric(metric), nodes(points.size()){
        for (int i = 0; i < (int)points.size(); i++) nodes[i].p = points[i], nodes[i].i = i;
        if (!nodes.empty()) build(0, nodes.size());
    }

    Neighbor nearest(const Point& q, int skip = -1) const {
        Neighbor best = {LLONG_MAX, -1};
        if (!nodes.empty()) search(0, q, skip, best);
        return best;
    }

    vector<int> report(long long x1, long long y1, long long x2, long long y2) const {
        vector<int> res;
        report(x1, y1, x2, y2, [&](int i){ res.push_back(i); });
        return res;
    }

    template <typename F>
    void report(long long x1, long long y1, long long x2, long long y2, F f) const {
        if (!nodes.empty()) collect(0, x1, y1, x2, y2, f);
    }

    int build(int l, int r){
        long long x1 = LLONG_MAX, x2 = LLONG_MIN, y1 = LLONG_MAX, y2 = LLONG_MIN;
        for (int k = l; k < r; k++){
            x1 = min(x1, nodes[k].p.x), x2 = max(x2, nodes[k].p.x);
            y1 = min(y1, nodes[k].p.y), y2 = max(y2, nodes[k].p.y);
        }

        int c = cells.size();
        cells.push_back({x1, x2, y1, y2, l, r, -1, -1});
        if (r - l <= LEAF || (x1 == x2 && y1 == y2)) return c;

        int m = split(l, r, x2 - x1 >= y2 - y1);
        int left = build(l, m);
        int right = build(m, r);
        cells[c].left = left, cells[c].right = right;
        return c;
    }

    /// One partition around the median of 15 samples is cheaper than nth_element, kept only for lopsided splits
    int split(int l, int r, bool by_x){
        auto key = [by_x](const Node& node){ return by_x ? node.p.x : node.p.y; };
        int len = r - l;
        long long sample[15];
        for (int j = 0; j < 15; j++) sample[j] = key(nodes[l + (2 * j + 1) * (len / 30)]);
        nth_element(sample, sample + 7, sample + 15);

        long long pivot = sample[7];
        auto first = nodes.begin() + l, last = nodes.begin() + r;
        auto mid = partition(first, last, [&](const Node& node){ return key(node) < pivot; });
        /// Points equal to the pivot join a too small left side, so duplicates share a cell that stops splitting
        if (mid - first <= len / 8) mid = partition(mid, last, [&](const Node& node){ return key(node) == pivot; });
        if (mid - first > len / 8 && last - mid > len / 8) return mid - nodes.begin();

        mid = first + len / 2;
        nth_element(first, mid, last, [&](const Node& a, const Node& b){ return key(a) < key(b); });
        return mid - nodes.begin();
    }

    void search(int c, const Point& q, int skip, Neighbor& best) const {
        const Cell& cell = cells[c];
        /// A cell of equal points is a leaf of any size, so its nearest point is found here without a scan
        if (cell.x1 == cell.x2 && cell.y1 == cell.y2){
            int k = nodes[cell.l].i == skip ? cell.l + 1 : cell.l;
            long long d = dist(cell.x1 - q.x, cell.y1 - q.y);
            if (k < cell.r && d < best.dist) best = {d, nodes[k].i};
            return;
        }

        if (cell.left < 0){
            scan(cell, q, skip, best);
            return;
        }

        long long dl = box_dist(cells[cell.left], q), dr = box_dist(cells[cell.right], q);
        if (dl <= dr){
            if (dl < best.dist) search(cell.left, q, skip, best);
            if (dr < best.dist) search(cell.right, q, skip, best);
        }
        else{
            if (dr < best.dist) search(cell.right, q, skip, best);
            if (dl < best.dist) search(cell.left, q, skip, best);
        }
    }

    /// The minimum first and its index only on improvement: the loop then has no branch on the distances
    void scan(const Cell& cell, const Point& q, int skip, Neighbor& best) const {
        long long leaf_best = LLONG_MAX;
        for (int k = cell.l; k < cell.r; k++){
            long long d = dist(nodes[k].p.x - q.x, nodes[k].p.y - q.y);
            leaf_best = min(leaf_best, nodes[k].i == skip ? LLONG_MAX : d);
        }
        if (leaf_best >= best.dist) return;

        int k = cell.l;
        while (nodes[k].i == skip || dist(nodes[k].p.x - q.x, nodes[k].p.y - q.y) != leaf_best) k++;
        best = {leaf_best, nodes[k].i};
    }

    long long box_dist(const Cell& cell, const Point& q) const {
        long long dx = max({cell.x1 - q.x, q.x - cell.x2, 0LL});
        long long dy = max({cell.y1 - q.y, q.y - cell.y2, 0LL});
        return dist(dx, dy);
    }

    long long dist(long long dx, long long dy) const {
        return metric == EUCLIDEAN ? dx * dx + dy * dy : abs(dx) + abs(dy);
    }

    template <typename F>
    void collect(int c, long long x1, long long y1, long long x2, long long y2, F& f) const {
        const Cell& cell = cells[c];
        if (cell.x2 < x1 || x2 < cell.x1 || cell.y2 < y1 || y2 < cell.y1) return;
        if (x1 <= cell.x1 && cell.x2 <= x2 && y1 <= cell.y1 && cell.y2 <= y2){
            for (int k = cell.l; k < cell.r; k++) f(nodes[k].i);
            return;
        }

        if (cell.left < 0){
            for (int k = cell.l; k < cell.r; k++){
                const Point& p = nodes[k].p;
                if (x1 <= p.x && p.x <= x2 && y1 <= p.y && p.y <= y2) f(nodes[k].i);
            }

            return;
        }

        collect(cell.left, x1, y1, x2, y2, f);
        collect(cell.right, x1, y1, x2, y2, f);
    }
};

int main(){
    /// The nearest point to the origin and the nearest other point of each input point
    vector<Point> points = {{3, 0}, {2, 2}, {7, 6}, {3, 0}};
    KDTree tree(points);
    assert(tree.nearest({0, 0}).dist == 8);
    assert(tree.nearest({0, 0}).i == 1);
    assert(tree.nearest(points[2], 2).dist == 41);
    assert(tree.nearest(points[2], 2).i == 1);
    assert(tree.nearest(points[0], 0).dist == 0);
    assert(tree.nearest(points[0], 0).i == 3);

    /// Every input point's nearest other point, visited in tree order
    vector<long long> nearest_other(points.size());
    for (const KDTree::Node& node : tree.nodes) nearest_other[node.i] = tree.nearest(node.p, node.i).dist;
    assert(nearest_other == vector<long long>({0, 5, 41, 0}));

    /// Manhattan distance prefers the point on the axis: |3| + |0| = 3 beats |2| + |2| = 4
    KDTree manhattan(points, KDTree::MANHATTAN);
    assert(manhattan.nearest({0, 0}).dist == 3);
    assert(manhattan.nearest(points[2], 2).dist == 9);
    assert(manhattan.nearest(points[2], 2).i == 1);

    /// Rectangles are inclusive on all four sides
    vector<int> inside = tree.report(2, 0, 3, 2);
    sort(inside.begin(), inside.end());
    assert(inside == vector<int>({0, 1, 3}));
    assert(tree.report(4, 0, 6, 5).empty());

    int count = 0;
    tree.report(2, 0, 3, 2, [&](int){ count++; });
    assert(count == 3);

    return 0;
}
