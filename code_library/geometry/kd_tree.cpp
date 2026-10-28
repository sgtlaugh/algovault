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
 *
 * Notes:
 *   - |x|, |y| <= 1e9 for points and queries, so a squared distance (at most 8e18) fits in a long long
 *   - nearest(points[i], i) is the nearest other point, at dist 0 when points[i] has a duplicate
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

    /// The subtree over positions [l, r) is rooted at (l + r) / 2, which holds the box of the whole subtree
    struct Node{
        Point p;
        int i;
        long long x1, x2, y1, y2;
    };

    Metric metric;
    vector<Node> nodes;

    KDTree(const vector<Point>& points, Metric metric = EUCLIDEAN) : metric(metric), nodes(points.size()){
        for (int i = 0; i < (int)points.size(); i++) nodes[i].p = points[i], nodes[i].i = i;
        build(0, nodes.size());
    }

    Neighbor nearest(const Point& q, int skip = -1) const {
        Neighbor best = {LLONG_MAX, -1};
        if (!nodes.empty()) search(0, nodes.size(), q, skip, best);
        return best;
    }

    vector<int> report(long long x1, long long y1, long long x2, long long y2) const {
        vector<int> res;
        collect(0, nodes.size(), x1, y1, x2, y2, res);
        return res;
    }

    void build(int l, int r){
        if (l >= r) return;

        long long x1 = LLONG_MAX, x2 = LLONG_MIN, y1 = LLONG_MAX, y2 = LLONG_MIN;
        for (int k = l; k < r; k++){
            x1 = min(x1, nodes[k].p.x), x2 = max(x2, nodes[k].p.x);
            y1 = min(y1, nodes[k].p.y), y2 = max(y2, nodes[k].p.y);
        }

        int m = (l + r) / 2;
        bool by_x = x2 - x1 >= y2 - y1;
        nth_element(nodes.begin() + l, nodes.begin() + m, nodes.begin() + r, [&](const Node& a, const Node& b){
            return by_x ? a.p.x < b.p.x : a.p.y < b.p.y;
        });
        nodes[m].x1 = x1, nodes[m].x2 = x2, nodes[m].y1 = y1, nodes[m].y2 = y2;

        build(l, m);
        build(m + 1, r);
    }

    void search(int l, int r, const Point& q, int skip, Neighbor& best) const {
        int m = (l + r) / 2;
        const Node& node = nodes[m];
        if (node.i != skip){
            long long d = dist(node.p.x - q.x, node.p.y - q.y);
            if (d < best.dist) best = {d, node.i};
        }

        long long dl = box_dist(l, m, q), dr = box_dist(m + 1, r, q);
        if (dl <= dr){
            if (dl < best.dist) search(l, m, q, skip, best);
            if (dr < best.dist) search(m + 1, r, q, skip, best);
        }
        else{
            if (dr < best.dist) search(m + 1, r, q, skip, best);
            if (dl < best.dist) search(l, m, q, skip, best);
        }
    }

    /// LLONG_MAX for an empty range, which no best distance can beat, so empty children are never visited
    long long box_dist(int l, int r, const Point& q) const {
        if (l >= r) return LLONG_MAX;

        const Node& node = nodes[(l + r) / 2];
        long long dx = max({node.x1 - q.x, q.x - node.x2, 0LL});
        long long dy = max({node.y1 - q.y, q.y - node.y2, 0LL});
        return dist(dx, dy);
    }

    long long dist(long long dx, long long dy) const {
        return metric == EUCLIDEAN ? dx * dx + dy * dy : abs(dx) + abs(dy);
    }

    void collect(int l, int r, long long x1, long long y1, long long x2, long long y2, vector<int>& res) const {
        if (l >= r) return;

        int m = (l + r) / 2;
        const Node& node = nodes[m];
        if (node.x2 < x1 || x2 < node.x1 || node.y2 < y1 || y2 < node.y1) return;
        if (x1 <= node.x1 && node.x2 <= x2 && y1 <= node.y1 && node.y2 <= y2){
            for (int k = l; k < r; k++) res.push_back(nodes[k].i);
            return;
        }

        if (x1 <= node.p.x && node.p.x <= x2 && y1 <= node.p.y && node.p.y <= y2) res.push_back(node.i);
        collect(l, m, x1, y1, x2, y2, res);
        collect(m + 1, r, x1, y1, x2, y2, res);
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

    return 0;
}
