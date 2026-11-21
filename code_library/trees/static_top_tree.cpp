/***
 *
 * Static Top Tree
 * A tree DP kept up to date under point updates, answered from the fixed root or from any vertex as the root
 *
 * Complexity: O(n) to build, O(log n) operation calls per update and per value
 *
 * StaticTopTree<Ops> tree(n, edges, ops, root = 0, roots = ANY_ROOT): vertices 0..n-1, edge e joins edges[e]
 * tree.update(v): recomputes after the data of vertex v changed
 * tree.update_edge(e): recomputes after the data of edge e changed
 * tree.value(): the DP of the whole tree rooted at root, a Path
 * tree.value(r): the DP of the whole tree rooted at r, a Path built as ops.add_vertex(..., r, -1) or vertex(r, -1)
 *
 * Notes:
 *   - A Path is a chain of vertices x_1..x_k, everything hanging off them, and the edge e joining x_1 to the outside
 *   - A Path's last vertex x_k may join the outside by an edge the Path does not hold, and e = -1 means none
 *   - A Point is everything hanging off one vertex through one edge, that edge included
 *   - Ops defines types Path and Point, both default constructible, and five calls, read from the data they hold:
 *       - vertex(v, e): the Path of v alone with edge e
 *       - add_vertex(x, v, e): the Path of v with edge e and the Point x hanging off v
 *       - add_edge(p): p as a Point hanging off the vertex its edge joins, p's last vertex joining nothing beyond
 *       - rake(a, b): two Points hanging off the same vertex as one, associative and commutative
 *       - compress(top, bottom): one Path, the edge of bottom joining its first vertex to the last vertex of top
 *   - value(r) joins Paths built from both ends of each chain, so the five calls must describe one consistent DP
 *   - ANY_ROOT keeps every Path in both directions; FIXED_ROOT keeps one, 15-30% faster, and value(r) then asserts
 *   - Ops is copied, so it should reach mutable data through pointers or references
 *   - Recursion depth is O(log n), deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Ops>
struct StaticTopTree{
    using Path = typename Ops::Path;
    using Point = typename Ops::Point;
    enum Kind{VERTEX, ADD_VERTEX, COMPRESS, ADD_EDGE, RAKE};
    enum Roots{ANY_ROOT, FIXED_ROOT};

    int n, root, top_node;
    Roots roots;
    Ops ops;
    vector<pair<int, int>> edges;
    vector<int> parent, parent_edge, heavy, heavy_edge, size, top, vertex_node;
    vector<char> kind;
    vector<int> left, right, up, vtx;
    vector<Path> fwd, rev;  /// rev is the same Path seen from its other end, edge heavy_edge instead of parent_edge
    vector<Point> point;

    StaticTopTree(int n, const vector<pair<int, int>>& edges, Ops ops, int root = 0, Roots roots = ANY_ROOT)
        : n(n), root(root), roots(roots), ops(ops), edges(edges), parent(n, -1), parent_edge(n, -1), heavy(n, -1), heavy_edge(n, -1),
          size(n, 1), top(n), vertex_node(n){
        assert(n >= 1 && (int)edges.size() == n - 1 && 0 <= root && root < n);

        vector<int> start(n + 1, 0), adj(2 * (n - 1));
        for (auto [u, v] : edges) start[u + 1]++, start[v + 1]++;
        for (int v = 0; v < n; v++) start[v + 1] += start[v];
        vector<int> pos(start.begin(), start.end() - 1);
        for (int e = 0; e < n - 1; e++) adj[pos[edges[e].first]++] = e, adj[pos[edges[e].second]++] = e;

        vector<int> order = {root};
        vector<char> seen(n, 0);
        seen[root] = 1;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int e = adj[j], v = edges[e].first ^ edges[e].second ^ u;
                if (seen[v]) continue;
                seen[v] = 1, parent[v] = u, parent_edge[v] = e;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i > 0; i--) size[parent[order[i]]] += size[order[i]];
        for (int i = 1; i < n; i++){
            int v = order[i], p = parent[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v, heavy_edge[p] = parent_edge[v];
        }

        vector<vector<int>> light(n);
        top[root] = root;
        for (int i = 1; i < n; i++){
            int v = order[i];
            top[v] = heavy[parent[v]] == v ? top[parent[v]] : v;
            if (heavy[parent[v]] != v) light[parent[v]].push_back(v);
        }
        top_node = build_chain(root, light);
        up[top_node] = -1;
    }

    void update(int v){
        for (int node = vertex_node[v]; node != -1; node = up[node]) recompute(node);
    }

    void update_edge(int e){
        update(edges[e].first), update(edges[e].second);
    }

    Path value() const{
        return fwd[top_node];
    }

    /// Climbs from r: each chain on the way gives the part below and above the vertex the climb entered it at
    Path value(int r){
        assert(roots == ANY_ROOT);

        struct Level{
            int v;
            bool has_light = false, has_below = false, has_above = false;
            Point light;
            Path below, above;  /// below is a forward Path, above a reversed one starting next to v
        };

        vector<Level> levels;
        Level cur;
        cur.v = r;
        int node = vertex_node[r];
        if (kind[node] == ADD_VERTEX) cur.has_light = true, cur.light = point[left[node]];

        for (int child = node, p = up[node]; p != -1; child = p, p = up[p]){
            if (kind[p] == COMPRESS && child == left[p]){
                cur.below = cur.has_below ? ops.compress(cur.below, fwd[right[p]]) : fwd[right[p]];
                cur.has_below = true;
            }
            else if (kind[p] == COMPRESS){
                cur.above = cur.has_above ? ops.compress(cur.above, rev[left[p]]) : rev[left[p]];
                cur.has_above = true;
            }
            else if (kind[p] == ADD_EDGE){
                levels.push_back(cur);
                cur = Level();
            }
            else if (kind[p] == RAKE){
                const Point& other = point[child == left[p] ? right[p] : left[p]];
                cur.light = cur.has_light ? ops.rake(cur.light, other) : other;
                cur.has_light = true;
            }
            else cur.v = vtx[p];
        }
        levels.push_back(cur);

        /// rest is everything beyond the top of the chain of levels[i], a Path whose edge joins that chain's top
        bool has_rest = levels.back().has_above;
        Path rest = levels.back().above;
        for (int i = (int)levels.size() - 2; i >= -1; i--){
            const Level& at = levels[i + 1];
            bool has_x = at.has_light;
            Point x = at.light;
            auto hang = [&](const Path& p){
                Point hanging = ops.add_edge(p);
                x = has_x ? ops.rake(x, hanging) : hanging;
                has_x = true;
            };
            if (at.has_below) hang(at.below);
            if (has_rest) hang(rest);

            if (i == -1) return has_x ? ops.add_vertex(x, at.v, -1) : ops.vertex(at.v, -1);

            int e = parent_edge[top[levels[i].v]];
            Path beyond = has_x ? ops.add_vertex(x, at.v, e) : ops.vertex(at.v, e);
            rest = levels[i].has_above ? ops.compress(levels[i].above, beyond) : beyond;
            has_rest = true;
        }

        return rest;
    }

    int build_chain(int head, const vector<vector<int>>& light){
        vector<int> nodes, weights;
        for (int v = head; v != -1; v = heavy[v]){
            vector<int> points, point_weights;
            for (int c : light[v]){
                points.push_back(add_node(ADD_EDGE, build_chain(c, light), -1, -1));
                point_weights.push_back(size[c]);
            }

            int node;
            if (points.empty()) node = add_node(VERTEX, -1, -1, v);
            else node = add_node(ADD_VERTEX, merge(points, point_weights, 0, points.size(), RAKE), -1, v);
            vertex_node[v] = node;
            nodes.push_back(node), weights.push_back(size[v] - (heavy[v] == -1 ? 0 : size[heavy[v]]));
        }

        return merge(nodes, weights, 0, nodes.size(), COMPRESS);
    }

    /// Splits at the weight midpoint, so a vertex lies under O(log n) nodes in all
    int merge(const vector<int>& nodes, const vector<int>& weights, int l, int r, Kind k){
        if (r - l == 1) return nodes[l];

        long long total = 0, prefix = 0;
        for (int i = l; i < r; i++) total += weights[i];
        int m = l + 1;
        for (prefix = weights[l]; m < r - 1 && 2 * (prefix + weights[m]) <= total; m++) prefix += weights[m];

        return add_node(k, merge(nodes, weights, l, m, k), merge(nodes, weights, m, r, k), -1);
    }

    int add_node(Kind k, int l, int r, int v){
        int node = kind.size();
        kind.push_back(k), left.push_back(l), right.push_back(r), up.push_back(-1), vtx.push_back(v);
        fwd.emplace_back(), rev.emplace_back(), point.emplace_back();
        if (l != -1) up[l] = node;
        if (r != -1) up[r] = node;
        recompute(node);
        return node;
    }

    void recompute(int node){
        int l = left[node], r = right[node], v = vtx[node];
        bool both = roots == ANY_ROOT;
        if (kind[node] == VERTEX){
            fwd[node] = ops.vertex(v, parent_edge[v]);
            if (both) rev[node] = ops.vertex(v, heavy_edge[v]);
        }
        else if (kind[node] == ADD_VERTEX){
            fwd[node] = ops.add_vertex(point[l], v, parent_edge[v]);
            if (both) rev[node] = ops.add_vertex(point[l], v, heavy_edge[v]);
        }
        else if (kind[node] == COMPRESS){
            fwd[node] = ops.compress(fwd[l], fwd[r]);
            if (both) rev[node] = ops.compress(rev[r], rev[l]);
        }
        else if (kind[node] == ADD_EDGE) point[node] = ops.add_edge(fwd[l]);
        else point[node] = ops.rake(point[l], point[r]);
    }
};

int main(){
    /***
     *             0
     *        3 /     \ 1
     *         1       2
     *    2 /    \ 7    \ 4
     *     3      4      5
     *  5 /
     *   6
    ***/
    vector<pair<int, int>> edges = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {3, 6}};
    vector<long long> weight = {3, 1, 2, 7, 4, 5};

    /// Distance to the farthest vertex: a Path keeps its length and its farthest vertex, both from its edge's far end
    struct Farthest{
        struct Path{ long long len = 0, far = 0; };
        struct Point{ long long far = 0; };
        vector<long long>* weight;

        long long w(int e) const{ return e == -1 ? 0 : (*weight)[e]; }
        Path vertex(int, int e) const{ return {w(e), w(e)}; }
        Path add_vertex(const Point& x, int, int e) const{ return {w(e), w(e) + x.far}; }
        Point add_edge(const Path& p) const{ return {p.far}; }
        Point rake(const Point& a, const Point& b) const{ return {max(a.far, b.far)}; }
        Path compress(const Path& top, const Path& bottom) const{
            return {top.len + bottom.len, max(top.far, top.len + bottom.far)};
        }
    };

    StaticTopTree<Farthest> tree(7, edges, Farthest{&weight});
    assert(tree.value().far == 10);   /// 0 -> 1 -> 3 -> 6
    assert(tree.value(4).far == 15);  /// 4 -> 1 -> 0 -> 2 -> 5
    assert(tree.value(2).far == 11);  /// 2 -> 0 -> 1 -> 3 -> 6

    /// Edge 3 joins 1 and 4: shortening it to 1 changes the farthest distance from 4
    weight[3] = 1;
    tree.update_edge(3);
    assert(tree.value().far == 10);
    assert(tree.value(4).far == 9);   /// 4 -> 1 -> 0 -> 2 -> 5

    /// Only ever rooted at 3, so each Path is kept in one direction
    StaticTopTree<Farthest> rooted(7, edges, Farthest{&weight}, 3, StaticTopTree<Farthest>::FIXED_ROOT);
    assert(rooted.value().far == 10);  /// 3 -> 1 -> 0 -> 2 -> 5

    return 0;
}
