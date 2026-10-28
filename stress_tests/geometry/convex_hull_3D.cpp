#include "../common.h"

#define main library_main
#include "../../code_library/geometry/convex_hull_3D.cpp"
#undef main

const long long C = 1000000000;

struct Vec{
    long long x, y, z;
    bool operator<(const Vec& o) const{ return tie(x, y, z) < tie(o.x, o.y, o.z); }
    bool operator==(const Vec& o) const{ return x == o.x && y == o.y && z == o.z; }
};

Vec normal(const Point3D& a, const Point3D& b, const Point3D& c){
    long long ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z, vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    return {uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx};
}

__int128 dot(const Vec& n, const Point3D& a, const Point3D& q){
    return (__int128)n.x * (q.x - a.x) + (__int128)n.y * (q.y - a.y) + (__int128)n.z * (q.z - a.z);
}

__int128 offset(const Vec& n, const Point3D& q){
    return (__int128)n.x * q.x + (__int128)n.y * q.y + (__int128)n.z * q.z;
}

Vec primitive(Vec n){
    long long g = gcd(gcd(llabs(n.x), llabs(n.y)), llabs(n.z));
    return {n.x / g, n.y / g, n.z / g};
}

struct Plane{
    Vec n;
    __int128 off;
    bool operator<(const Plane& o) const{ return tie(n, off) < tie(o.n, o.off); }
};

Plane plane_of(const vector<Point3D>& pts, int a, int b, int c){
    Vec n = primitive(normal(pts[a], pts[b], pts[c]));
    return {n, offset(n, pts[a])};
}

long long cross2(pair<long long, long long> o, pair<long long, long long> a, pair<long long, long long> b){
    return (a.first - o.first) * (b.second - o.second) - (a.second - o.second) * (b.first - o.first);
}

long long hull_area2(vector<pair<long long, long long>> p){
    sort(p.begin(), p.end());
    p.erase(unique(p.begin(), p.end()), p.end());
    if (p.size() < 3) return 0;

    vector<pair<long long, long long>> h(2 * p.size());
    int k = 0;
    for (int i = 0; i < (int)p.size(); i++){
        while (k >= 2 && cross2(h[k - 2], h[k - 1], p[i]) <= 0) k--;
        h[k++] = p[i];
    }
    for (int i = (int)p.size() - 2, lower = k + 1; i >= 0; i--){
        while (k >= lower && cross2(h[k - 2], h[k - 1], p[i]) <= 0) k--;
        h[k++] = p[i];
    }

    long long res = 0;
    for (int i = 1; i + 1 < k; i++) res += cross2(h[0], h[i], h[i + 1]);
    return res;
}

/// Every face has positive area and no point strictly above it, and the faces form a closed oriented surface
void check_surface(const vector<Point3D>& pts, const ConvexHull3D& hull, int support_faces){
    int f = hull.faces.size();
    map<pair<int, int>, int> edges;
    set<int> verts;
    for (auto& t : hull.faces){
        for (int i = 0; i < 3; i++){
            assert(0 <= t[i] && t[i] < (int)pts.size());
            verts.insert(t[i]);
            edges[{t[i], t[(i + 1) % 3]}]++;
        }
        Vec n = normal(pts[t[0]], pts[t[1]], pts[t[2]]);
        assert(n.x || n.y || n.z);
    }
    for (auto& [e, cnt] : edges){
        assert(cnt == 1);
        assert(edges.count({e.second, e.first}));
    }
    assert((int)verts.size() - (int)edges.size() / 2 + f == 2);

    for (int k = 0; k < min(f, support_faces); k++){
        auto& t = hull.faces[support_faces >= f ? k : stress::rand_int(0, f - 1)];
        Vec n = normal(pts[t[0]], pts[t[1]], pts[t[2]]);
        for (const Point3D& q : pts) assert(dot(n, pts[t[0]], q) <= 0);
    }
}

__int128 volume6(const vector<Point3D>& pts, const ConvexHull3D& hull){
    __int128 res = 0;
    for (auto& t : hull.faces) res -= dot(normal(pts[t[0]], pts[t[1]], pts[t[2]]), pts[t[0]], pts[hull.faces[0][0]]);
    return res;
}

bool close(long double got, long double want){
    return fabsl(got - want) <= 1e-12L * max((long double)1, fabsl(want));
}

/// A triple spans a facet iff no point is above it; checks facet planes, the area each covers, volume and area
void check_brute(const vector<Point3D>& pts){
    int n = pts.size();
    ConvexHull3D hull(pts);

    set<Plane> facets;
    bool solid = false;
    for (int a = 0; a < n; a++){
        for (int b = 0; b < n; b++){
            for (int c = 0; c < n; c++){
                Vec nm = normal(pts[a], pts[b], pts[c]);
                if (!nm.x && !nm.y && !nm.z) continue;

                bool support = true, flat = true;
                for (int q = 0; q < n && support; q++){
                    __int128 s = dot(nm, pts[a], pts[q]);
                    support = s <= 0;
                    flat &= s == 0;
                }
                solid |= !flat;
                if (support && !flat) facets.insert(plane_of(pts, a, b, c));
            }
        }
    }

    if (!solid){
        assert(hull.faces.empty());
        assert(hull.volume() == 0);
        assert(hull.area() == 0);
        return;
    }
    check_surface(pts, hull, INT_MAX);

    map<Plane, Vec> covered;
    for (auto& t : hull.faces){
        Plane pl = plane_of(pts, t[0], t[1], t[2]);
        assert(facets.count(pl));
        Vec nm = normal(pts[t[0]], pts[t[1]], pts[t[2]]);
        Vec& sum = covered[pl];
        sum = {sum.x + nm.x, sum.y + nm.y, sum.z + nm.z};
    }
    assert(covered.size() == facets.size());

    __int128 vol6 = 0;
    long double area = 0;
    for (const Plane& pl : facets){
        Vec nm = pl.n;
        long long ax = llabs(nm.x), ay = llabs(nm.y), az = llabs(nm.z);
        int axis = ax >= ay && ax >= az ? 0 : ay >= az ? 1 : 2;
        vector<pair<long long, long long>> proj;
        for (int q = 0; q < n; q++){
            if (offset(nm, pts[q]) != pl.off) continue;
            if (axis == 0) proj.push_back({pts[q].y, pts[q].z});
            if (axis == 1) proj.push_back({pts[q].z, pts[q].x});
            if (axis == 2) proj.push_back({pts[q].x, pts[q].y});
        }

        long long a2 = hull_area2(proj), big = max({ax, ay, az});
        assert(a2 % big == 0);
        long long t = a2 / big;
        assert((covered[pl] == Vec{t * nm.x, t * nm.y, t * nm.z}));
        vol6 += t * pl.off;
        area += t * sqrtl((long double)ax * ax + (long double)ay * ay + (long double)az * az) / 2;
    }

    assert(volume6(pts, hull) == vol6);
    assert(close(hull.volume(), (long double)vol6 / 6));
    assert(close(hull.area(), area));
}

vector<Point3D> random_points(int n, long long lo, long long hi){
    vector<Point3D> pts;
    for (int i = 0; i < n; i++) pts.push_back({stress::rand_int(lo, hi), stress::rand_int(lo, hi), stress::rand_int(lo, hi)});
    return pts;
}

/// Points on k random planes z = a x + b y + c, so most facets hold many coplanar points
vector<Point3D> plane_points(int n, int k, long long range){
    vector<array<long long, 3>> planes;
    for (int i = 0; i < k; i++) planes.push_back({stress::rand_int(-2, 2), stress::rand_int(-2, 2), stress::rand_int(-3, 3)});
    vector<Point3D> pts;
    for (int i = 0; i < n; i++){
        auto& p = planes[stress::rand_int(0, k - 1)];
        long long x = stress::rand_int(-range, range), y = stress::rand_int(-range, range);
        int perm = stress::rand_int(0, 2);
        Point3D q = {x, y, p[0] * x + p[1] * y + p[2]};
        if (perm == 1) q = {q.z, q.x, q.y};
        if (perm == 2) q = {q.y, q.z, q.x};
        pts.push_back(q);
    }

    return pts;
}

/// Directions rounded onto a sphere of radius r: nearly every point is a vertex, the worst case for hull size
vector<Point3D> sphere_points(int n, long long r){
    vector<Point3D> pts;
    while ((int)pts.size() < n){
        long long x = stress::rand_int(-1e6, 1e6), y = stress::rand_int(-1e6, 1e6), z = stress::rand_int(-1e6, 1e6);
        long long len2 = x * x + y * y + z * z;
        if (len2 < 1e10 || len2 > 1e12) continue;

        long double len = sqrtl(len2);
        pts.push_back({llroundl(x / len * r), llroundl(y / len * r), llroundl(z / len * r)});
    }

    return pts;
}

/// Every lattice point of a cube: each side is a flat facet holding (side + 1)^2 coplanar points
void check_box(int side){
    vector<Point3D> pts;
    for (int x = 0; x <= side; x++){
        for (int y = 0; y <= side; y++){
            for (int z = 0; z <= side; z++) pts.push_back({x, y, z});
        }
    }
    shuffle(pts.begin(), pts.end(), stress::rng());

    ConvexHull3D hull(pts);
    check_surface(pts, hull, INT_MAX);
    assert(volume6(pts, hull) == 6LL * side * side * side);
    assert(hull.area() == 6LL * side * side);
}

/// The corners of the coordinate bound, points on one side and inside: 6 V = 6 (2e9)^3 overflows long long
void check_bound(int n){
    vector<Point3D> pts = random_points(n, -C, C);
    for (int m = 0; m < 8; m++) pts.push_back({m & 1 ? C : -C, m & 2 ? C : -C, m & 4 ? C : -C});
    for (int i = 0; i < n; i++) pts.push_back({stress::rand_int(0, 1) ? C : -C, stress::rand_int(-C, C), stress::rand_int(-C, C)});
    shuffle(pts.begin(), pts.end(), stress::rng());

    ConvexHull3D hull(pts);
    check_surface(pts, hull, 50);
    __int128 side = 2 * C;
    assert(volume6(pts, hull) == 6 * side * side * side);
    assert(close(hull.volume(), (long double)side * side * side));
    assert(close(hull.area(), 6 * (long double)side * side));
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(0, it % 10 ? 9 : 16);
        int kind = it % 6;
        vector<Point3D> pts;
        if (kind == 0) pts = random_points(n, -1, 1);
        if (kind == 1) pts = random_points(n, 0, 2);
        if (kind == 2) pts = random_points(n, -C, C);
        if (kind == 3) pts = plane_points(n, stress::rand_int(1, 3), 3);
        if (kind == 4) pts = random_points(n, C - 2, C);
        if (kind == 5){
            pts = random_points(n, -5, 5);
            for (auto& p : pts) p.z = it % 3 == 0 ? p.x + 2 * p.y : p.x * 3;
            for (auto& p : pts) p.y = it % 3 == 2 ? -2 * p.x : p.y;
        }
        check_brute(pts);
    }

    for (int it = 0; it < 6; it++){
        vector<Point3D> pts = it % 3 == 0 ? random_points(2000, -1000, 1000) : it % 3 == 1 ? sphere_points(2000, 300) : plane_points(2000, 4, 30);
        ConvexHull3D hull(pts);
        check_surface(pts, hull, INT_MAX);
    }
    check_box(12);
    check_bound(2000);

    /// Sorted by height, so a hull that inserts in input order instead of shuffling takes minutes, not a second
    vector<Point3D> big = sphere_points(30000, C);
    sort(big.begin(), big.end(), [](const Point3D& p, const Point3D& q){ return p.z < q.z; });
    ConvexHull3D hull(big);
    check_surface(big, hull, 100);

    return 0;
}
