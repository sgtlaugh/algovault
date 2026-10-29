/***
 *
 * 3D Convex Hull
 * Triangulated convex hull of integer points in 3D, with its volume and surface area
 *
 * Complexity: expected O(n log n) time, memory peaks near 150 MB for 1e5 points that are all vertices
 *
 * ConvexHull3D hull(points): the hull of points, faces empty when every point lies on one plane
 * hull.area(): surface area
 * hull.volume(): volume
 * hull.faces: triangles {a, b, c} as indices into points, counter-clockwise seen from outside
 *
 * Notes:
 *   - |x|, |y|, |z| <= 1e9, every orientation test is exact in __int128 (64-bit GCC or Clang)
 *   - No input point lies strictly outside the plane of any face, and no face has zero area
 *   - Duplicate points are allowed, a vertex is any one of its copies
 *   - A flat facet comes as coplanar triangles, maybe with vertices inside it or on its edges; count facets by plane
 *   - The insertion order is random, so the triangulation of a flat facet can change between runs
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point3D{
    long long x, y, z;
};

struct ConvexHull3D{
    struct Face{
        array<int, 3> v, next;  /// next[i] is the face across the edge v[i] -> v[(i + 1) % 3]
        long long nx, ny, nz;
        vector<int> outside;  /// unadded points strictly above this face, emptied when the face dies
        bool dead = false;
    };

    vector<Point3D> points;
    vector<array<int, 3>> faces;

    ConvexHull3D(const vector<Point3D>& points) : points(points){
        int n = points.size(), a = 0, b = 1, c = 0, d = 0;
        auto same = [&](int i){
            return points[i].x == points[a].x && points[i].y == points[a].y && points[i].z == points[a].z;
        };
        auto collinear = [&](int i){
            Face f = make_face(a, b, i);
            return !f.nx && !f.ny && !f.nz;
        };
        while (b < n && same(b)) b++;
        while (b < n && c < n && collinear(c)) c++;
        if (b >= n || c >= n) return;

        Face abc = make_face(a, b, c);
        while (d < n && side(abc, d) == 0) d++;
        if (d == n) return;

        if (side(abc, d) > 0) swap(b, c);
        vector<Face> all;
        int base[4][3] = {{a, b, c}, {a, d, b}, {b, d, c}, {c, d, a}};
        int across[4][3] = {{1, 2, 3}, {3, 2, 0}, {1, 3, 0}, {2, 1, 0}};
        for (int f = 0; f < 4; f++){
            all.push_back(make_face(base[f][0], base[f][1], base[f][2]));
            all[f].next = {across[f][0], across[f][1], across[f][2]};
        }

        vector<vector<int>> sees(n);
        vector<int> order;
        for (int p = 0; p < n; p++){
            if (p == a || p == b || p == c || p == d) continue;
            order.push_back(p);
            for (int f = 0; f < 4; f++){
                if (side(all[f], p) > 0) all[f].outside.push_back(p), sees[p].push_back(f);
            }
        }
        /// Random insertion order is what makes the expected total conflict work O(n log n)
        static mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        shuffle(order.begin(), order.end(), rng);

        vector<int> start(n), stamp(n, -1), visible, fresh;
        for (int p : order){
            visible.clear();
            for (int f : sees[p]){
                if (!all[f].dead) all[f].dead = true, visible.push_back(f);
            }
            sees[p] = vector<int>();
            if (visible.empty()) continue;

            fresh.clear();
            for (int f : visible){
                for (int i = 0; i < 3; i++){
                    int g = all[f].next[i];
                    if (all[g].dead) continue;

                    int u = all[f].v[i], w = all[f].v[(i + 1) % 3], h = all.size();
                    all.push_back(make_face(u, w, p));
                    all[h].next[0] = g;
                    for (int j = 0; j < 3; j++){
                        if (all[g].v[j] == w) all[g].next[j] = h;
                    }
                    start[u] = h;
                    fresh.push_back(h);

                    for (const vector<int>* list : {&all[f].outside, &all[g].outside}){
                        for (int q : *list){
                            if (stamp[q] == h || side(all[h], q) <= 0) continue;

                            stamp[q] = h;
                            all[h].outside.push_back(q);
                            sees[q].push_back(h);
                        }
                    }
                }
            }

            for (int h : fresh) all[h].next[1] = start[all[h].v[1]];
            for (int h : fresh) all[all[h].next[1]].next[2] = h;
            for (int f : visible) all[f].outside = vector<int>();
        }

        for (const Face& f : all){
            if (!f.dead) faces.push_back(f.v);
        }
    }

    long double area() const{
        long double res = 0;
        for (auto& f : faces){
            Face g = make_face(f[0], f[1], f[2]);
            res += sqrtl((long double)g.nx * g.nx + (long double)g.ny * g.ny + (long double)g.nz * g.nz);
        }

        return res / 2;
    }

    long double volume() const{
        __int128 res = 0;
        for (auto& f : faces) res -= side(make_face(f[0], f[1], f[2]), faces[0][0]);

        return (long double)res / 6;
    }

    Face make_face(int a, int b, int c) const{
        const Point3D &p = points[a], &q = points[b], &r = points[c];
        long long ux = q.x - p.x, uy = q.y - p.y, uz = q.z - p.z, vx = r.x - p.x, vy = r.y - p.y, vz = r.z - p.z;
        Face f;
        f.v = {a, b, c};
        f.nx = uy * vz - uz * vy, f.ny = uz * vx - ux * vz, f.nz = ux * vy - uy * vx;
        return f;
    }

    /// Six times the signed volume of the face's triangle and q, positive when q is outside
    __int128 side(const Face& f, int q) const{
        const Point3D &p = points[f.v[0]], &r = points[q];
        return (__int128)f.nx * (r.x - p.x) + (__int128)f.ny * (r.y - p.y) + (__int128)f.nz * (r.z - p.z);
    }
};

int main(){
    /// A cube of side 2 with its center and a repeated corner: two triangles per square side, 8 vertices
    vector<Point3D> cube = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {2, 2, 0}, {0, 0, 2}, {2, 0, 2}, {0, 2, 2}, {2, 2, 2}};
    cube.push_back({1, 1, 1});
    cube.push_back({2, 2, 2});
    ConvexHull3D hull(cube);
    assert(hull.faces.size() == 12);
    assert(hull.volume() == 8);
    assert(hull.area() == 24);

    ConvexHull3D flat({{0, 0, 5}, {3, 0, 5}, {0, 4, 5}, {1, 1, 5}});
    assert(flat.faces.empty());
    assert(flat.volume() == 0);

    return 0;
}
