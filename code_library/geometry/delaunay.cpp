/***
 *
 * Delaunay Triangulation
 * Triangulates points so that no point lies strictly inside any triangle's circumcircle
 *
 * Complexity: O(n log n)
 *
 * Delaunay dt(points): the Delaunay triangulation of points in any order, duplicates allowed
 * dt.triangles: index triples {a, b, c} into points, counter-clockwise with a the smallest
 * dt.edges: index pairs {u, v}, every edge of the triangulation once
 *
 * Notes:
 *   - |x|, |y| <= 1e9, requires __int128
 *   - Collinear points give no triangles, and edges join consecutive points along the line
 *   - Four or more cocircular points get any one of their valid triangulations
 *   - Duplicates: triangles and edges use the copy with the smallest index, edges also join every other copy u to it
 *   - edges contains a Euclidean MST: Kruskal from minimum_spanning_tree.cpp with weight dist2(u, v) finds it
 *   - Defines its own Point: pasted next to another geometry file, merge the Point blocks first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;
};

/// Guibas-Stolfi divide and conquer on quad-edges, adapted from KACTL FastDelaunay.h (MIT licence)
struct Delaunay{
    vector<Point> pts;
    vector<int> id;
    /// Directed edge e runs org[e] -> org[e ^ 1], onext and oprev link the edges leaving org[e] in ccw order
    vector<int> org, onext, oprev, free_edges;
    vector<array<int, 3>> triangles;
    vector<array<int, 2>> edges;

    Delaunay(const vector<Point>& points){
        int n = points.size();
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){
            return make_tuple(points[a].x, points[a].y, a) < make_tuple(points[b].x, points[b].y, b);
        });

        for (int i : order){
            const Point& p = points[i];
            if (!pts.empty() && pts.back().x == p.x && pts.back().y == p.y) edges.push_back({i, id.back()});
            else pts.push_back(p), id.push_back(i);
        }

        if (pts.size() >= 2) build(0, pts.size());
        collect();
    }

    long long cross(int o, int a, int b) const{
        return (pts[a].x - pts[o].x) * (pts[b].y - pts[o].y) - (pts[a].y - pts[o].y) * (pts[b].x - pts[o].x);
    }

    /// True if d is strictly inside the circle through a, b, c in ccw order
    bool in_circle(int a, int b, int c, int d) const{
        auto lift = [&](int u){
            long long dx = pts[u].x - pts[d].x, dy = pts[u].y - pts[d].y;
            return (__int128)(dx * dx + dy * dy);
        };
        return lift(a) * cross(d, b, c) + lift(b) * cross(d, c, a) + lift(c) * cross(d, a, b) > 0;
    }

    int dest(int e) const{
        return org[e ^ 1];
    }

    int lnext(int e) const{
        return oprev[e ^ 1];
    }

    /// Returns {ccw hull edge out of the leftmost point, cw hull edge out of the rightmost point} of pts[l, r)
    pair<int, int> build(int l, int r){
        if (r - l == 2){
            int a = make_edge(l, l + 1);
            return {a, a ^ 1};
        }
        if (r - l == 3){
            int a = make_edge(l, l + 1), b = make_edge(l + 1, l + 2);
            splice(a ^ 1, b);
            long long turn = cross(l, l + 1, l + 2);
            if (turn == 0) return {a, b ^ 1};
            int c = connect(b, a);
            if (turn > 0) return {a, b ^ 1};
            return {c ^ 1, c};
        }

        auto [ldo, ldi] = build(l, (l + r) / 2);
        auto [rdi, rdo] = build((l + r) / 2, r);

        while (true){
            if (cross(org[ldi], dest(ldi), org[rdi]) > 0) ldi = lnext(ldi);
            else if (cross(org[rdi], dest(rdi), org[ldi]) < 0) rdi = onext[rdi ^ 1];
            else break;
        }

        int base = connect(rdi ^ 1, ldi);
        if (org[ldi] == org[ldo]) ldo = base ^ 1;
        if (org[rdi] == org[rdo]) rdo = base;

        /// Zips the two halves bottom up, base runs from the right half to the left half along the lower boundary
        auto valid = [&](int e){ return cross(dest(e), dest(base), org[base]) > 0; };
        while (true){
            int lcand = onext[base ^ 1];
            while (valid(lcand) && in_circle(dest(base), org[base], dest(lcand), dest(onext[lcand]))){
                int next = onext[lcand];
                remove_edge(lcand);
                lcand = next;
            }

            int rcand = oprev[base];
            while (valid(rcand) && in_circle(dest(base), org[base], dest(rcand), dest(oprev[rcand]))){
                int next = oprev[rcand];
                remove_edge(rcand);
                rcand = next;
            }

            if (!valid(lcand) && !valid(rcand)) break;
            if (!valid(lcand) || (valid(rcand) && in_circle(dest(lcand), org[lcand], org[rcand], dest(rcand)))){
                base = connect(rcand, base ^ 1);
            }
            else base = connect(base ^ 1, lcand ^ 1);
        }

        return {ldo, rdo};
    }

    int make_edge(int u, int v){
        int e = org.size();
        if (!free_edges.empty()){
            e = free_edges.back();
            free_edges.pop_back();
        }
        else{
            org.resize(e + 2);
            onext.resize(e + 2);
            oprev.resize(e + 2);
        }

        org[e] = u, org[e ^ 1] = v;
        onext[e] = oprev[e] = e;
        onext[e ^ 1] = oprev[e ^ 1] = e ^ 1;
        return e;
    }

    /// Merges the rings around org[a] and org[b] if they differ, else splits them, as in Guibas-Stolfi
    void splice(int a, int b){
        int an = onext[a], bn = onext[b];
        onext[a] = bn, onext[b] = an;
        oprev[bn] = a, oprev[an] = b;
    }

    /// Adds the edge dest(a) -> org[b], with a's left face to its left and b's left face to its right
    int connect(int a, int b){
        int e = make_edge(dest(a), org[b]);
        splice(e, lnext(a));
        splice(e ^ 1, b);
        return e;
    }

    void remove_edge(int e){
        splice(e, oprev[e]);
        splice(e ^ 1, oprev[e ^ 1]);
        free_edges.push_back(e & ~1);
    }

    void collect(){
        int m = org.size();
        vector<bool> alive(m, true);
        for (int e : free_edges) alive[e] = alive[e ^ 1] = false;

        for (int e = 0; e < m; e += 2){
            if (alive[e]) edges.push_back({id[org[e]], id[org[e ^ 1]]});
        }

        /// Each triangle is a ccw face of three edges, reported once from its smallest edge; the outer face turns cw
        for (int a = 0; a < m; a++){
            if (!alive[a]) continue;
            int b = lnext(a), c = lnext(b);
            if (lnext(c) == a && a < b && a < c && cross(org[a], org[b], org[c]) > 0){
                array<int, 3> t = {id[org[a]], id[org[b]], id[org[c]]};
                rotate(t.begin(), min_element(t.begin(), t.end()), t.end());
                triangles.push_back(t);
            }
        }
    }
};

int main(){
    /// A square with its center: four triangles around the center, eight edges
    vector<Point> points = {{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 2}};
    Delaunay dt(points);
    assert(dt.triangles.size() == 4);
    assert(dt.edges.size() == 8);
    for (auto [a, b, c] : dt.triangles) assert(a == 4 || b == 4 || c == 4);

    Delaunay tri({{0, 0}, {5, 1}, {1, 3}});
    assert((tri.triangles == vector<array<int, 3>>{{0, 1, 2}}));

    Delaunay line({{6, 6}, {0, 0}, {3, 3}});
    assert(line.triangles.empty());
    assert(line.edges.size() == 2);

    Delaunay twins({{1, 2}, {7, 7}, {1, 2}});
    assert(twins.edges.size() == 2);
    assert(find(twins.edges.begin(), twins.edges.end(), array<int, 2>{2, 0}) != twins.edges.end());

    return 0;
}
