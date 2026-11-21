#include "../common.h"

#define main library_main
#include "../../code_library/geometry/geometry_3D.cpp"
#undef main

/// Brute forces minimise over the shapes' parameters by nested ternary search, valid since each distance is convex
template <typename F>
double ternary_min(double lo, double hi, F f, int steps = 45){
    for (int it = 0; it < steps; it++){
        double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        if (f(m1) <= f(m2)) hi = m2;
        else lo = m1;
    }

    return f((lo + hi) / 2);
}

Point3 at(const Point3& a, const Point3& b, double s){
    return a + (b - a) * s;
}

Point3 at(const Point3& a, const Point3& b, const Point3& c, double u, double v){
    return a + (b - a) * u + (c - a) * v;
}

double brute_point_segment(const Point3& p, const Point3& a, const Point3& b){
    return ternary_min(0, 1, [&](double s){ return dist(p, at(a, b, s)); });
}

double brute_segment_segment(const Point3& a, const Point3& b, const Point3& c, const Point3& d){
    return ternary_min(0, 1, [&](double s){ return brute_point_segment(at(a, b, s), c, d); });
}

double brute_point_triangle(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    return ternary_min(0, 1, [&](double u){
        return ternary_min(0, 1 - u, [&](double v){ return dist(p, at(a, b, c, u, v)); });
    });
}

double brute_segment_triangle(const Point3& p, const Point3& q, const Point3& a, const Point3& b, const Point3& c){
    return ternary_min(0, 1, [&](double s){ return brute_point_triangle(at(p, q, s), a, b, c); });
}

long long rnd(int lim){
    return stress::rand_int(-lim, lim);
}

Point3 rand_point(int lim){
    return {(double)rnd(lim), (double)rnd(lim), (double)rnd(lim)};
}

/// Orientation of p against plane abc, exact for small integer coordinates
long long side(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    long long ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z, vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    long long wx = p.x - a.x, wy = p.y - a.y, wz = p.z - a.z;
    return wx * (uy * vz - uz * vy) + wy * (uz * vx - ux * vz) + wz * (ux * vy - uy * vx);
}

void check_close(double got, double want, const char* what){
    if (fabs(got - want) > 1e-6){
        fprintf(stderr, "%s: got %.12f, want %.12f\n", what, got, want);
        abort();
    }
}

/// Picks points that make the degenerate cases common: repeats, collinear triples and parallel segments
Point3 shape_point(const vector<Point3>& made, int lim){
    int kind = stress::rand_int(0, 9);
    if (kind == 0 && !made.empty()) return made[stress::rand_int(0, made.size() - 1)];
    if (kind == 1 && made.size() >= 2){
        Point3 a = made[stress::rand_int(0, made.size() - 1)], b = made[stress::rand_int(0, made.size() - 1)];
        return a + (b - a) * (double)stress::rand_int(-2, 3);
    }

    return rand_point(lim);
}

void distances(int iterations, int lim){
    for (int it = 0; it < iterations; it++){
        vector<Point3> pts;
        for (int i = 0; i < 6; i++) pts.push_back(shape_point(pts, lim));
        if (stress::rand_int(0, 4) == 0) pts[3] = pts[2] + (pts[1] - pts[0]) * (double)stress::rand_int(-2, 2);

        Point3 p = pts[0], q = pts[1], a = pts[2], b = pts[3], c = pts[4], x = pts[5];
        check_close(dist_point_segment(x, a, b), brute_point_segment(x, a, b), "point segment");
        check_close(dist_segment_segment(p, q, a, b), brute_segment_segment(p, q, a, b), "segment segment");
        check_close(dist_point_triangle(x, a, b, c), brute_point_triangle(x, a, b, c), "point triangle");
        check_close(dist_segment_triangle(p, q, a, b, c), brute_segment_triangle(p, q, a, b, c), "segment triangle");

        if (dist(a, b) > 0){
            double line = ternary_min(-1000, 1000, [&](double s){ return dist(x, at(a, b, s)); }, 90);
            check_close(dist_point_line(x, a, b), line, "point line");
        }

        long long s = side(x, a, b, c);
        if (norm(cross(b - a, c - a)) > 0){
            double plane = dist_point_plane(x, a, b, c);
            assert((plane > 1e-9) == (s > 0));
            assert((plane < -1e-9) == (s < 0));
            double far = ternary_min(-500, 500, [&](double u){
                return ternary_min(-500, 500, [&](double v){ return dist(x, at(a, b, c, u, v)); }, 90);
            }, 90);
            check_close(fabs(plane), far, "point plane");
        }
    }
}

/// A segment through the inside of a triangle, or ending on its boundary, is at distance exactly 0
void touching(int iterations){
    for (int it = 0; it < iterations; it++){
        Point3 a = rand_point(20), b = rand_point(20), c = rand_point(20);
        if (norm(cross(b - a, c - a)) == 0) continue;
        double u = stress::rand_int(0, 100) / 100.0, v = stress::rand_int(0, 100) / 100.0;
        if (u + v > 1) u = 1 - u, v = 1 - v;
        Point3 hit = at(a, b, c, u, v), dir = rand_point(20);
        if (norm(dir) == 0) continue;
        assert(dist_segment_triangle(hit - dir, hit + dir, a, b, c) < 1e-9);
        assert(dist_segment_triangle(hit, hit + dir, a, b, c) < 1e-9);
        assert(dist_point_triangle(hit, a, b, c) < 1e-9);
    }
}

void plane_intersections(int iterations){
    for (int it = 0; it < iterations; it++){
        Point3 a = rand_point(30), b = rand_point(30), c = rand_point(30), p = rand_point(30), q = rand_point(30);
        if (norm(cross(b - a, c - a)) == 0 || dist(p, q) == 0) continue;

        long long sp = side(p, a, b, c), sq = side(q, a, b, c);
        pair<int, Point3> res = line_plane_intersection(p, q, a, b, c);
        if (sp == sq){
            assert(res.first == (sp == 0 ? -1 : 0));
            continue;
        }

        assert(res.first == 1);
        double scale = max({norm(p), norm(q), norm(a), norm(b), norm(c), 1.0});
        assert(fabs(dist_point_plane(res.second, a, b, c)) < 1e-9 * scale);
        assert(dist_point_line(res.second, p, q) < 1e-9 * scale);
    }

    /// Parallel lines and lines in the plane of z = 5
    Point3 a{0, 0, 5}, b{7, 0, 5}, c{0, 3, 5};
    assert(line_plane_intersection({1, 1, 4}, {9, -3, 4}, a, b, c).first == 0);
    assert(line_plane_intersection({1, 1, 5}, {9, -3, 5}, a, b, c).first == -1);
}

/// Boxes, tetrahedra and triangular prisms with known volume, shifted far away and with faces listed in random order
void volumes(int iterations){
    for (int it = 0; it < iterations; it++){
        Point3 o = rand_point(1000);
        double w = stress::rand_int(1, 50), h = stress::rand_int(1, 50), d = stress::rand_int(1, 50);
        vector<Point3> box;
        for (int i = 0; i < 8; i++) box.push_back(o + Point3{i & 1 ? w : 0, i & 2 ? h : 0, i & 4 ? d : 0});
        vector<vector<int>> faces = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
        shuffle(faces.begin(), faces.end(), stress::rng());
        for (auto& f : faces) rotate(f.begin(), f.begin() + stress::rand_int(0, 3), f.end());
        check_close(polyhedron_volume(box, faces) / (w * h * d), 1, "box");

        vector<Point3> tet = {rand_point(30), rand_point(30), rand_point(30), rand_point(30)};
        long long s6 = side(tet[3], tet[0], tet[1], tet[2]);
        if (s6 == 0) continue;
        vector<vector<int>> tf = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {0, 2, 3}};
        if (s6 > 0) for (auto& f : tf) swap(f[1], f[2]);
        check_close(polyhedron_volume(tet, tf), llabs(s6) / 6.0, "tetrahedron");

        /// Prism: a triangle at z = 0 and its copy at z = H, volume = base area * H
        Point3 t0{(double)rnd(20), (double)rnd(20), 0}, t1{(double)rnd(20), (double)rnd(20), 0};
        Point3 t2{(double)rnd(20), (double)rnd(20), 0};
        double base = cross(t1 - t0, t2 - t0).z / 2, H = stress::rand_int(1, 20);
        if (base == 0) continue;
        if (base < 0) swap(t1, t2), base = -base;
        vector<Point3> pr = {t0, t1, t2, t0 + Point3{0, 0, H}, t1 + Point3{0, 0, H}, t2 + Point3{0, 0, H}};
        vector<vector<int>> pf = {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
        check_close(polyhedron_volume(pr, pf), base * H, "prism");
    }
}

/// Long, almost parallel segments, where solving the two-line system without the parallel guard loses about 1e-6
void near_parallel(int iterations){
    for (int it = 0; it < iterations; it++){
        Point3 dir{1e6 * stress::rand_int(1, 1000) / 1000.0 + 1, (double)rnd(6), (double)rnd(4)};
        Point3 a{0, 0, 0}, off{(double)rnd(4), (double)rnd(4), 1};
        Point3 c = off + dir * (stress::rand_int(0, 99) / 100.0);
        double e = stress::rand_int(0, 999) * 1e-9;
        Point3 b = a + dir, d = c + dir * 0.7 + Point3{0, e, e};
        if (it & 1){
            /// A long thin X: the lines cross mid-way 1e-3 apart at an angle of about 1e-5, far from every endpoint
            double k = stress::rand_int(1, 9);
            c = {0, k, 1e-3}, d = {dir.x, -k, 1e-3}, b = {dir.x, 0, 0};
        }

        auto gap = [&](long double s, long double t){
            long double x = a.x + (b.x - a.x) * s - c.x - (d.x - c.x) * t;
            long double y = a.y + (b.y - a.y) * s - c.y - (d.y - c.y) * t;
            long double z = a.z + (b.z - a.z) * s - c.z - (d.z - c.z) * t;
            return sqrtl(x * x + y * y + z * z);
        };
        long double want = ternary_min(0, 1, [&](double s){
            return (double)ternary_min(0, 1, [&](double t){ return (double)gap(s, t); }, 100);
        }, 100);
        assert(fabs(dist_segment_segment(a, b, c, d) - want) < 1e-7);
    }
}

int main(){
    distances(stress::scaled(30), 2);
    distances(stress::scaled(30), 10);
    near_parallel(stress::scaled(300));
    touching(stress::scaled(3000));
    plane_intersections(stress::scaled(5000));
    volumes(stress::scaled(2000));
    return 0;
}
