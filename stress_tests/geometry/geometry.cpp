#include "../common.h"

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

Point random_point(long long range){
    return {stress::rand_int(-range, range), stress::rand_int(-range, range)};
}

/// Lattice points of a segment are a + k * step for step = (b - a) / g, a point is on it iff it is one of them
bool brute_on_segment(Point p, Point a, Point b){
    long long dx = b.x - a.x, dy = b.y - a.y, g = __gcd(llabs(dx), llabs(dy));
    if (g == 0) return p == a;
    for (long long k = 0; k <= g; k++){
        if (a.x + k * (dx / g) == p.x && a.y + k * (dy / g) == p.y) return true;
    }
    return false;
}

/// Intersection by solving the two line equations exactly, overlap of collinear segments by projection
bool brute_intersect(Point a, Point b, Point c, Point d){
    if (a == b) return brute_on_segment(a, c, d);
    if (c == d) return brute_on_segment(c, a, b);

    __int128 r_x = b.x - a.x, r_y = b.y - a.y, s_x = d.x - c.x, s_y = d.y - c.y;
    __int128 denom = r_x * s_y - r_y * s_x, qp_x = c.x - a.x, qp_y = c.y - a.y;
    if (denom == 0){
        if (qp_x * r_y - qp_y * r_x != 0) return false;
        __int128 rr = r_x * r_x + r_y * r_y;
        __int128 t0 = qp_x * r_x + qp_y * r_y, t1 = t0 + s_x * r_x + s_y * r_y;
        return max(min(t0, t1), (__int128)0) <= min(max(t0, t1), rr);
    }

    __int128 t = qp_x * s_y - qp_y * s_x, u = qp_x * r_y - qp_y * r_x;
    if (denom < 0) denom = -denom, t = -t, u = -u;
    return 0 <= t && t <= denom && 0 <= u && u <= denom;
}

/// Winding number from the sum of signed angles, exact boundary test by lattice enumeration
int brute_point_in_polygon(const vector<Point>& poly, Point p){
    int n = poly.size();
    for (int i = 0; i < n; i++){
        if (brute_on_segment(p, poly[i], poly[(i + 1) % n])) return 0;
    }

    long double total = 0;
    for (int i = 0; i < n; i++){
        long double ax = poly[i].x - p.x, ay = poly[i].y - p.y, bx = poly[(i + 1) % n].x - p.x, by = poly[(i + 1) % n].y - p.y;
        total += atan2l(ax * by - ay * bx, ax * bx + ay * by);
    }
    return fabsl(total) > 3 ? 1 : -1;
}

/// Simple: every point used once, non-adjacent edges disjoint, adjacent edges meet only at their shared vertex
bool is_simple(const vector<Point>& pts, const vector<int>& order){
    int n = pts.size();
    vector<int> sorted = order;
    sort(sorted.begin(), sorted.end());
    for (int i = 0; i < n; i++){
        if (sorted[i] != i) return false;
    }

    for (int i = 0; i < n; i++){
        Point a = pts[order[i]], b = pts[order[(i + 1) % n]];
        for (int j = i + 1; j < n; j++){
            Point c = pts[order[j]], d = pts[order[(j + 1) % n]];
            bool adjacent = j == i + 1 || (i == 0 && j == n - 1);
            if (!adjacent){
                if (brute_intersect(a, b, c, d)) return false;
            }
            else if (j == i + 1){
                if (brute_on_segment(a, c, d) || brute_on_segment(d, a, b)) return false;
            }
            else if (brute_on_segment(b, c, d) || brute_on_segment(c, a, b)) return false;
        }
    }
    return true;
}

long double brute_dist(Point p, Point a, Point b){
    long double lo = 0, hi = 1;
    auto f = [&](long double t){ return hypotl(a.x + t * (b.x - a.x) - p.x, a.y + t * (b.y - a.y) - p.y); };
    for (int it = 0; it < 100; it++){
        long double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        if (f(m1) < f(m2)) hi = m2;
        else lo = m1;
    }
    return min({f(lo), f(0), f(1)});
}

vector<Point> brute_convex_hull(vector<Point> pts){
    vector<Point> hull;
    for (Point p : pts){
        bool vertex = true;
        for (Point q : pts){
            for (Point r : pts){
                if (q == p || r == p || q == r) continue;
                if (cross(q, r, p) == 0 && on_segment(p, q, r)) vertex = false;
            }
        }
        if (vertex) hull.push_back(p);
    }

    vector<Point> res;
    for (Point p : hull){
        bool extreme = false;
        for (Point q : hull){
            if (q == p) continue;
            bool all_left = true;
            for (Point r : pts){
                if (cross(p, q, r) < 0) all_left = false;
            }
            if (all_left) extreme = true;
        }
        if (extreme) res.push_back(p);
    }

    if (res.size() < 3) return {};
    Point o = *min_element(res.begin(), res.end(), [](Point u, Point v){ return make_pair(u.y, u.x) < make_pair(v.y, v.x); });
    sort(res.begin(), res.end(), [&](Point u, Point v){
        if (u == o) return !(v == o);
        if (v == o) return false;
        long long c = cross(o, u, v);
        return c != 0 ? c > 0 : dot(o, u, u) < dot(o, v, v);
    });
    return res;
}

PointF to_f(Point p){
    return {(double)p.x, (double)p.y};
}

bool close_to(long double got, long double expected){
    return fabsl(got - expected) <= 1e-9L * (1 + fabsl(expected));
}

/// p = s + t * (e - s) with 0 <= t <= 1, decided exactly in __int128
bool exact_on(Point p, Point s, Point e){
    __int128 rx = e.x - s.x, ry = e.y - s.y, qx = p.x - s.x, qy = p.y - s.y;
    if (rx == 0 && ry == 0) return qx == 0 && qy == 0;
    __int128 t = qx * rx + qy * ry;
    return qx * ry - qy * rx == 0 && 0 <= t && t <= rx * rx + ry * ry;
}

using Points = vector<pair<long double, long double>>;

/// The common points of two segments from the exact line equations: a fraction for a crossing, the overlap interval otherwise
Points brute_segment_points(Point a, Point b, Point c, Point d){
    auto at = [&](__int128 num, __int128 den){ return make_pair(a.x + (long double)((b.x - a.x) * num) / den, a.y + (long double)((b.y - a.y) * num) / den); };
    if (a == b && c == d) return a == c ? Points{{a.x, a.y}} : Points{};
    if (a == b) return exact_on(a, c, d) ? Points{{a.x, a.y}} : Points{};
    if (c == d) return exact_on(c, a, b) ? Points{{c.x, c.y}} : Points{};

    __int128 rx = b.x - a.x, ry = b.y - a.y, sx = d.x - c.x, sy = d.y - c.y, qx = c.x - a.x, qy = c.y - a.y;
    __int128 den = rx * sy - ry * sx;
    if (den != 0){
        if (!brute_intersect(a, b, c, d)) return {};
        return {at(qx * sy - qy * sx, den)};
    }
    if (qx * ry - qy * rx != 0) return {};

    __int128 rr = rx * rx + ry * ry, t0 = qx * rx + qy * ry, t1 = t0 + sx * rx + sy * ry;
    __int128 lo = max(min(t0, t1), (__int128)0), hi = min(max(t0, t1), rr);
    if (lo > hi) return {};
    if (lo == hi) return {at(lo, rr)};
    return {at(lo, rr), at(hi, rr)};
}

bool same_point_sets(const vector<PointF>& got, const Points& expected){
    if (got.size() != expected.size()) return false;
    auto has = [&](long double x, long double y){
        for (const PointF& p : got){
            if (close_to(p.x, x) && close_to(p.y, y)) return true;
        }
        return false;
    };
    for (auto [x, y] : expected){
        if (!has(x, y)) return false;
    }
    return true;
}

/// Lines, segments, projection, reflection and distance against exact fractions, orientation against the __int128 sign
void check_lines(Point a, Point b, Point c, Point d, Point p){
    PointF fa = to_f(a), fb = to_f(b), fc = to_f(c), fd = to_f(d), fp = to_f(p);
    auto sign = [](__int128 v){ return (v > 0) - (v < 0); };
    assert(orientation(fa, fb, fc) == sign((__int128)(b.x - a.x) * (c.y - a.y) - (__int128)(b.y - a.y) * (c.x - a.x)));
    assert(same_point_sets(segment_intersection(fa, fb, fc, fd), brute_segment_points(a, b, c, d)));
    if (a == b || c == d) return;

    __int128 rx = b.x - a.x, ry = b.y - a.y, sx = d.x - c.x, sy = d.y - c.y, qx = c.x - a.x, qy = c.y - a.y;
    __int128 den = rx * sy - ry * sx;
    auto [kind, point] = line_intersection(fa, fb, fc, fd);
    if (den == 0) assert(kind == (qx * ry - qy * rx == 0 ? -1 : 0));
    else{
        __int128 t = qx * sy - qy * sx;
        assert(kind == 1);
        assert(close_to(point.x, (long double)(a.x * den + rx * t) / den) && close_to(point.y, (long double)(a.y * den + ry * t) / den));
    }

    __int128 rr = rx * rx + ry * ry, t = (__int128)(p.x - a.x) * rx + (__int128)(p.y - a.y) * ry;
    long double px = (long double)(a.x * rr + rx * t) / rr, py = (long double)(a.y * rr + ry * t) / rr;
    PointF proj = project(fp, fa, fb), refl = reflect(fp, fa, fb);
    assert(close_to(proj.x, px) && close_to(proj.y, py));
    assert(close_to(refl.x, 2 * px - p.x) && close_to(refl.y, 2 * py - p.y));

    long double dist = hypotl(p.x - px, p.y - py);
    int side = sign(rx * (p.y - a.y) - ry * (p.x - a.x));
    assert(close_to(dist_point_line(fp, fa, fb), side * dist) && (side == 0 || (dist_point_line(fp, fa, fb) > 0) == (side > 0)));
}

/// Same direction by reduced vectors, otherwise atan2l in [0, 2 pi), with the exact __int128 cross only for angles closer than 1e-12
bool brute_angle_less(Point a, Point b){
    bool za = a.x == 0 && a.y == 0, zb = b.x == 0 && b.y == 0;
    if (za || zb) return za && !zb;

    long long ga = __gcd(llabs(a.x), llabs(a.y)), gb = __gcd(llabs(b.x), llabs(b.y));
    if (a.x / ga == b.x / gb && a.y / ga == b.y / gb) return false;

    long double pi = acosl(-1.0L), ta = atan2l(a.y, a.x), tb = atan2l(b.y, b.x);
    if (ta < 0) ta += 2 * pi;
    if (tb < 0) tb += 2 * pi;
    if (fabsl(ta - tb) > 1e-12) return ta < tb;
    return (__int128)a.x * b.y - (__int128)a.y * b.x > 0;
}

/// Winding number of a point that is never on the boundary, by the sum of signed angles
bool brute_inside(const vector<PointF>& poly, long double x, long double y){
    long double total = 0;
    for (int i = 0, n = poly.size(); i < n; i++){
        long double ax = poly[i].x - x, ay = poly[i].y - y, bx = poly[(i + 1) % n].x - x, by = poly[(i + 1) % n].y - y;
        total += atan2l(ax * by - ay * bx, ax * bx + ay * by);
    }
    return fabsl(total) > 3;
}

long double area_f(const vector<PointF>& poly){
    long double s = 0;
    for (int i = 0, n = poly.size(); i < n; i++) s += (long double)poly[i].x * poly[(i + 1) % n].y - (long double)poly[(i + 1) % n].x * poly[i].y;
    return s / 2;
}

/// Monte Carlo area of poly cut by the half-plane, the complementary cut, and for convex polygons point membership
void check_cut(const vector<Point>& poly, Point a, Point b, bool convex){
    vector<PointF> fpoly;
    for (Point p : poly) fpoly.push_back(to_f(p));
    PointF fa = to_f(a), fb = to_f(b);
    vector<PointF> kept = cut_polygon(fpoly, fa, fb), rest = cut_polygon(fpoly, fb, fa);
    bool any_left = false, any_right = false;
    for (Point p : poly) any_left |= cross(a, b, p) > 0, any_right |= cross(a, b, p) < 0;
    assert(kept.empty() == !any_left && rest.empty() == !any_right);
    long double total = area_f(fpoly), scale = 1 + fabsl(total);
    assert(fabsl(area_f(kept) + area_f(rest) - total) < 1e-9 * scale);

    long long lx = LLONG_MAX, hx = LLONG_MIN, ly = LLONG_MAX, hy = LLONG_MIN;
    for (Point p : poly) lx = min(lx, p.x), hx = max(hx, p.x), ly = min(ly, p.y), hy = max(hy, p.y);
    const int samples = 1500;
    int hits = 0;
    for (int i = 0; i < samples; i++){
        long double x = lx + (hx - lx) * (stress::rand_int(0, 1LL << 40) / (long double)(1LL << 40));
        long double y = ly + (hy - ly) * (stress::rand_int(0, 1LL << 40) / (long double)(1LL << 40));
        bool expected = brute_inside(fpoly, x, y) && (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x) > 0;
        hits += expected;
        if (!convex) continue;

        bool near_edge = fabsl((b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x)) < 1e-6 * hypotl(b.x - a.x, b.y - a.y);
        for (int j = 0, n = poly.size(); j < n; j++){
            long double ux = poly[j].x, uy = poly[j].y, vx = poly[(j + 1) % n].x - ux, vy = poly[(j + 1) % n].y - uy;
            long double t = max(0.0L, min(1.0L, ((x - ux) * vx + (y - uy) * vy) / (vx * vx + vy * vy)));
            near_edge |= hypotl(ux + t * vx - x, uy + t * vy - y) < 1e-6;
        }
        if (near_edge) continue;
        bool got = kept.size() >= 3;
        for (int j = 0, m = kept.size(); j < m && got; j++){
            const PointF& u = kept[j];
            const PointF& v = kept[(j + 1) % m];
            got = (v.x - u.x) * (y - u.y) - (v.y - u.y) * (x - u.x) > -1e-9;
        }
        assert(got == expected);
    }

    /// The spread comes from the claimed area, the observed hit rate understates it when only a few samples land
    long double box = (long double)(hx - lx) * (hy - ly), claimed = min(1.0L, max(0.0L, area_f(kept) / box));
    long double sd = box * sqrtl(max(claimed * (1 - claimed), 1.0L / samples) / samples);
    assert(fabsl(area_f(kept) - box * hits / samples) <= 6 * sd);

    if (convex){
        for (int j = 0, m = kept.size(); j < m; j++){
            assert(!(kept[j].x == kept[(j + 1) % m].x && kept[j].y == kept[(j + 1) % m].y));
            assert(cross(kept[j], kept[(j + 1) % m], kept[(j + 2) % m]) >= -1e-9);
        }
    }
}

/// Centroid from Green's theorem on the trapezoids under each edge, exact in __int128
pair<long double, long double> brute_centroid(const vector<Point>& poly){
    __int128 a6 = 0, mx6 = 0, my6 = 0;
    for (int i = 0, n = poly.size(); i < n; i++){
        __int128 x1 = poly[i].x, y1 = poly[i].y, dx = poly[(i + 1) % n].x - x1, dy = poly[(i + 1) % n].y - y1;
        a6 -= 3 * dx * (2 * y1 + dy);
        mx6 -= dx * (6 * x1 * y1 + 3 * (x1 * dy + y1 * dx) + 2 * dx * dy);
        my6 -= dx * (3 * y1 * y1 + 3 * y1 * dy + dy * dy);
    }
    return {(long double)mx6 / a6, (long double)my6 / a6};
}

int main(){
    const long long E9 = 1000000000;

    for (long long it = 0; it < stress::scaled(200000); it++){
        long long range = it % 3 ? 6 : 30;
        Point a = random_point(range), b = random_point(range), c = random_point(range), d = random_point(range), p = random_point(range);
        if (it % 5 == 0) b = a;
        assert(on_segment(p, a, b) == brute_on_segment(p, a, b));
        assert(segments_intersect(a, b, c, d) == brute_intersect(a, b, c, d));
        assert(abs(dist_point_segment(p, a, b) - brute_dist(p, a, b)) < 1e-6);
    }

    /// Coordinates drawn from a few values near +-1e9 so touching and collinear cases still occur at full range
    auto near_bound = [&](){
        switch (stress::rand_int(0, 3)){
            case 0: return -E9 + stress::rand_int(0, 2);
            case 1: return E9 - stress::rand_int(0, 2);
            case 2: return 0LL;
            default: return stress::rand_int(-E9, E9);
        }
    };

    for (long long it = 0; it < stress::scaled(20000); it++){
        Point a{near_bound(), near_bound()}, b{near_bound(), near_bound()}, c{near_bound(), near_bound()}, d{near_bound(), near_bound()};
        Point p{near_bound(), near_bound()};
        if (a == b || c == d) continue;  /// brute_on_segment walks every lattice point, too slow at this range
        assert(segments_intersect(a, b, c, d) == brute_intersect(a, b, c, d));
        long double expected = brute_dist(p, a, b);
        assert(fabsl(dist_point_segment(p, a, b) - expected) < 1e-6 + 1e-12 * expected);
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        Point a = random_point(E9), b = random_point(E9), p = random_point(E9);
        long long k = stress::rand_int(0, 7);
        Point on = {a.x + (b.x - a.x) / 7 * k, a.y + (b.y - a.y) / 7 * k};
        bool exact = (b.x - a.x) % 7 == 0 && (b.y - a.y) % 7 == 0;
        if (exact) assert(on_segment(on, a, b));
        assert(on_segment(a, a, b) && on_segment(b, a, b));
        assert(cross(a, b, p) == -cross(b, a, p));
    }

    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(3, 12);
        long long range = it % 2 ? 5 : 40;
        vector<Point> pts;
        set<pair<long long, long long>> used;
        while ((int)pts.size() < n){
            Point p = random_point(range);
            if (used.insert({p.x, p.y}).second) pts.push_back(p);
        }

        bool collinear = true;
        for (int i = 2; i < n; i++) collinear &= cross(pts[0], pts[1], pts[i]) == 0;
        if (collinear) continue;

        vector<int> order = simple_polygon(pts);
        assert(is_simple(pts, order));
        vector<Point> poly;
        for (int i : order) poly.push_back(pts[i]);

        __int128 shoelace = 0;
        for (int i = 0; i < n; i++) shoelace += (__int128)poly[i].x * poly[(i + 1) % n].y - (__int128)poly[(i + 1) % n].x * poly[i].y;
        assert(area2(poly) == shoelace);

        for (int q = 0; q < 30; q++){
            Point p = random_point(range + 2);
            assert(point_in_polygon(poly, p) == brute_point_in_polygon(poly, p));
            vector<Point> rev(poly.rbegin(), poly.rend());
            assert(point_in_polygon(rev, p) == brute_point_in_polygon(poly, p));
        }

        vector<Point> hull = brute_convex_hull(pts);
        if (hull.size() >= 3){
            for (int q = 0; q < 30; q++){
                Point p = random_point(range + 2);
                assert(point_in_convex_polygon(hull, p) == brute_point_in_polygon(hull, p));
            }
            for (Point v : hull) assert(point_in_convex_polygon(hull, v) == 0);
        }
    }

    /// Convex polygons with extra lattice vertices on their edges, rotated so vertex 0 can sit in the middle of an edge
    long long flat_starts = 0;
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(3, 10);
        long long range = stress::rand_int(2, 8);
        set<pair<long long, long long>> distinct;
        for (int i = 0; i < n; i++){
            Point p = random_point(range);
            distinct.insert({p.x, p.y});
        }
        vector<Point> pts;
        for (auto [x, y] : distinct) pts.push_back({x, y});
        vector<Point> hull = brute_convex_hull(pts);
        if (hull.size() < 3 || area2(hull) == 0) continue;

        vector<Point> poly;
        for (int i = 0; i < (int)hull.size(); i++){
            Point a = hull[i], b = hull[(i + 1) % hull.size()];
            long long g = __gcd(llabs(b.x - a.x), llabs(b.y - a.y));
            for (long long k = 0; k < g; k++){
                if (k == 0 || stress::rand_int(0, 1)) poly.push_back({a.x + k * ((b.x - a.x) / g), a.y + k * ((b.y - a.y) / g)});
            }
        }
        rotate(poly.begin(), poly.begin() + stress::rand_int(0, poly.size() - 1), poly.end());

        int m = poly.size();
        flat_starts += cross(poly[m - 1], poly[0], poly[1]) == 0;
        for (long long x = -range - 1; x <= range + 1; x++){
            for (long long y = -range - 1; y <= range + 1; y++){
                int got = point_in_convex_polygon(poly, {x, y});
                assert(got == brute_point_in_polygon(hull, {x, y}) && got == point_in_polygon(poly, {x, y}));
            }
        }
    }

    assert(flat_starts > 0);

    /// A spiral whose fan triangulation from vertex 0 has partial sums beyond long long, the total area fits
    vector<Point> spiral = {
        {-909090900, -999999990}, {999999990, -999999990}, {999999990, 999999990}, {-999999990, 999999990},
        {-999999990, -636363630}, {636363630, -636363630}, {636363630, 636363630}, {-636363630, 636363630},
        {-636363630, -272727270}, {272727270, -272727270}, {272727270, 181818180}, {90909090, 181818180},
        {90909090, -90909090}, {-454545450, -90909090}, {-454545450, 454545450}, {454545450, 454545450},
        {454545450, -454545450}, {-818181810, -454545450}, {-818181810, 818181810}, {818181810, 818181810},
        {818181810, -818181810}, {-909090900, -818181810},
    };
    __int128 spiral_area = 0, partial = 0, peak = 0;
    for (int i = 0; i < (int)spiral.size(); i++){
        Point u = spiral[i], v = spiral[(i + 1) % spiral.size()];
        spiral_area += (__int128)u.x * v.y - (__int128)v.x * u.y;
    }
    for (int i = 1; i + 1 < (int)spiral.size(); i++) partial += cross(spiral[0], spiral[i], spiral[i + 1]), peak = max(peak, partial < 0 ? -partial : partial);
    assert(peak > LLONG_MAX && area2(spiral) == spiral_area);

    /// Coordinates at the 1e9 boundary
    for (long long it = 0; it < stress::scaled(2000); it++){
        vector<Point> box = {{-E9, -E9}, {E9, -E9}, {E9, E9}, {-E9, E9}};
        Point p = random_point(E9);
        assert(point_in_convex_polygon(box, p) == point_in_polygon(box, p));
        assert(point_in_polygon(box, p) == ((llabs(p.x) == E9 || llabs(p.y) == E9) ? 0 : 1));
        assert(area2(box) == 8 * E9 * E9);
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        PointF c{stress::rand_int(-1000, 1000) / 7.0, stress::rand_int(-1000, 1000) / 7.0};
        PointF p{stress::rand_int(-1000, 1000) / 3.0, stress::rand_int(-1000, 1000) / 3.0};
        double deg = stress::rand_int(-720000, 720000) / 1000.0;
        complex<double> z = complex<double>(p.x - c.x, p.y - c.y) * polar(1.0, deg * acos(-1.0) / 180.0);
        PointF r = rotate(c, p, deg);
        assert(abs(r.x - (c.x + z.real())) < 1e-7 && abs(r.y - (c.y + z.imag())) < 1e-7);

        PointF u{p.x - c.x, p.y - c.y}, v{r.x - c.x, r.y - c.y};
        if (hypot(u.x, u.y) > 1e-3){
            double expected = fmod(-deg, 360.0);
            if (expected < 0) expected += 360.0;
            double got = clockwise_angle(u, v), diff = fabs(got - expected);
            assert(min(diff, 360.0 - diff) < 1e-6 && 0 <= got && got < 360.0);
        }
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        PointF u{(double)stress::rand_int(-1000000, 1000000), (double)stress::rand_int(-1000000, 1000000)};
        if (u.x == 0 && u.y == 0) continue;
        PointF v = rotate({0, 0}, u, stress::rand_int(0, 1) ? 1e-12 : 0.0);
        double got = clockwise_angle(u, v);
        assert(0 <= got && got < 360.0);
    }

    assert(clockwise_angle({1, 0}, {1, 1e-300}) < 360.0 && clockwise_angle({1, 0}, {1, -1e-300}) < 360.0);

    for (long long it = 0; it < stress::scaled(20000); it++){
        long double lat1 = stress::rand_int(-90000000, 90000000) / 1e6L, lon1 = stress::rand_int(-180000000, 180000000) / 1e6L;
        long double lat2 = it % 4 ? stress::rand_int(-90000000, 90000000) / 1e6L : lat1 + stress::rand_int(-10, 10) / 1e9L;
        long double lon2 = it % 4 ? stress::rand_int(-180000000, 180000000) / 1e6L : lon1 + stress::rand_int(-10, 10) / 1e9L;
        long double k = acosl(-1.0L) / 180;
        auto xyz = [&](long double lat, long double lon){ return array<long double, 3>{cosl(lat * k) * cosl(lon * k), cosl(lat * k) * sinl(lon * k), sinl(lat * k)}; };
        auto p = xyz(lat1, lon1), q = xyz(lat2, lon2);
        long double chord = sqrtl(powl(p[0] - q[0], 2) + powl(p[1] - q[1], 2) + powl(p[2] - q[2], 2));
        long double expected = 2 * asinl(min(1.0L, chord / 2));
        long double got = great_circle_distance(lat1, lon1, lat2, lon2, 1);
        assert(fabsl(got - expected) < 1e-9);
    }

    /// The haversine term of about 4% of antipodal pairs rounds to 1 + 1 ulp, the distance must still be pi * radius
    for (long long it = 0; it < stress::scaled(20000); it++){
        long double lat = stress::rand_int(-90000000, 90000000) / 1e6L, lon = stress::rand_int(-180000000, 0) / 1e6L;
        long double got = great_circle_distance(lat, lon, -lat, lon + 180, 6371);
        assert(fabsl(got - 6371 * acosl(-1.0L)) < 1e-9 * 6371);
    }

    /// Near-collinear floating inputs, which the integer-derived cases never produce, must fall inside the EPS tolerance
    assert(orientation({0, 0}, {1, 0}, {2, 1e-12}) == 0 && orientation({0, 0}, {1e6, 0}, {2e6, 1e-4}) == 0);
    assert(orientation({0, 0}, {1, 0}, {2, 1e-6}) == 1 && orientation({0, 0}, {1, 0}, {2, -1e-6}) == -1);
    assert(line_intersection({0, 0}, {1, 0}, {0, 1}, {1, 1 + 1e-12}).first == 0);

    /// Small coordinates hit every degenerate case: shared endpoints, collinear overlaps, zero-length segments
    for (long long it = 0; it < stress::scaled(150000); it++){
        long long range = it % 3 ? 4 : 25;
        Point a = random_point(range), b = random_point(range), c = random_point(range), d = random_point(range), p = random_point(range);
        if (it % 7 == 0) b = a;
        if (it % 11 == 0) d = c;
        check_lines(a, b, c, d, p);
    }

    /// The 1e4 boundary where integer input must still be decided exactly: long segments, points on or one unit off the line a b
    const long long E4 = 10000;
    auto edge_value = [&](){
        switch (stress::rand_int(0, 2)){
            case 0: return E4 - stress::rand_int(0, 1);
            case 1: return -E4 + stress::rand_int(0, 1);
            default: return stress::rand_int(-E4, E4);
        }
    };
    auto near_line = [&](Point a, Point b){
        long long dx = b.x - a.x, dy = b.y - a.y, g = __gcd(llabs(dx), llabs(dy));
        for (;;){
            long long k = stress::rand_int(-g / 2, g + g / 2);
            Point q{a.x + k * (dx / g), a.y + k * (dy / g)};
            if (stress::rand_int(0, 2) == 0) (stress::rand_int(0, 1) ? q.x : q.y) += stress::rand_int(0, 1) ? 1 : -1;
            if (llabs(q.x) <= E4 && llabs(q.y) <= E4) return q;
        }
    };
    long long collinear_hits = 0, near_misses = 0;
    for (long long it = 0; it < stress::scaled(60000); it++){
        Point a{edge_value(), edge_value()}, b{edge_value(), edge_value()};
        if (a == b) continue;
        Point c = near_line(a, b), d = it % 2 ? near_line(a, b) : Point{edge_value(), edge_value()}, p = near_line(a, b);
        long long side = cross(a, b, c);
        collinear_hits += side == 0, near_misses += llabs(side) == llabs(b.x - a.x) || llabs(side) == llabs(b.y - a.y);
        check_lines(a, b, c, d, p);
        check_lines(c, d, a, b, p);
    }

    assert(collinear_hits > 1000 && near_misses > 1000);

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long range = it % 2 ? 3 : 60;
        vector<Point> v(stress::rand_int(1, 30));
        for (Point& q : v) q = random_point(range);
        for (Point u : v){
            for (Point w : v) assert(angle_less(u, w) == brute_angle_less(u, w));
        }
        sort(v.begin(), v.end(), angle_less);
        for (int i = 0; i + 1 < (int)v.size(); i++) assert(!brute_angle_less(v[i + 1], v[i]));
    }

    /// Components at the 1e9 boundary, where a.x * b.y alone is close to 1e18
    auto big = [&](){
        switch (stress::rand_int(0, 5)){
            case 0: return E9 - stress::rand_int(0, 2);
            case 1: return -E9 + stress::rand_int(0, 2);
            case 2: return stress::rand_int(-2, 2);
            default: return stress::rand_int(-E9, E9);
        }
    };
    for (long long it = 0; it < stress::scaled(200000); it++){
        Point u{big(), big()}, w{big(), big()};
        assert(angle_less(u, w) == brute_angle_less(u, w));
    }

    /// Vectors p - o for p, o within 1e9 reach 2e9 per component, where the cross product is near 8e18
    auto big2 = [&](){
        switch (stress::rand_int(0, 5)){
            case 0: return 2 * E9 - stress::rand_int(0, 2);
            case 1: return -2 * E9 + stress::rand_int(0, 2);
            case 2: return stress::rand_int(-2, 2);
            default: return stress::rand_int(-2 * E9, 2 * E9);
        }
    };
    long long near_limit = 0;
    for (long long it = 0; it < stress::scaled(200000); it++){
        Point u{big2(), big2()}, w{big2(), big2()};
        assert(angle_less(u, w) == brute_angle_less(u, w));
        Point o{near_bound(), near_bound()}, p{near_bound(), near_bound()}, q{near_bound(), near_bound()};
        assert(angle_less(p - o, q - o) == brute_angle_less(p - o, q - o));
        near_limit += (__int128)u.x * w.y - (__int128)u.y * w.x > (__int128)7 * E9 * E9;
    }

    assert(near_limit > 0);

    /// Combs with three or more teeth tips on the line and the rest strictly on one side, the kept part on that side is empty
    for (long long it = 0; it < stress::scaled(2000); it++){
        int teeth = stress::rand_int(2, 6);
        long long ox = stress::rand_int(-20, 20), oy = stress::rand_int(-20, 20);
        vector<Point> comb = {{ox, oy - 10}, {ox + 2 * teeth, oy - 10}};
        for (int i = teeth; i >= 0; i--){
            comb.push_back({ox + 2 * i, oy});
            if (i > 0) comb.push_back({ox + 2 * i - 1, oy - stress::rand_int(1, 9)});
        }
        Point a{ox + 2 * stress::rand_int(0, teeth), oy}, b{ox + 2 * stress::rand_int(0, teeth), oy};
        if (!(a == b)) check_cut(comb, a, b, false);
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(3, 10);
        long long range = it % 2 ? 6 : 20;
        vector<Point> pts;
        set<pair<long long, long long>> used;
        while ((int)pts.size() < n){
            Point p = random_point(range);
            if (used.insert({p.x, p.y}).second) pts.push_back(p);
        }

        /// Every fourth line runs through two vertices, so vertices and whole edges lie exactly on it
        auto pick_line = [&](const vector<Point>& poly){
            if (it % 4 == 0) return make_pair(poly[stress::rand_int(0, poly.size() - 1)], poly[stress::rand_int(0, poly.size() - 1)]);
            return make_pair(random_point(range + 3), random_point(range + 3));
        };
        vector<Point> hull = brute_convex_hull(pts);
        if (hull.size() >= 3){
            auto [a, b] = pick_line(hull);
            if (!(a == b)) check_cut(hull, a, b, true);
        }

        bool collinear = true;
        for (int i = 2; i < n; i++) collinear &= cross(pts[0], pts[1], pts[i]) == 0;
        if (collinear) continue;
        vector<Point> poly;
        for (int i : simple_polygon(pts)) poly.push_back(pts[i]);
        auto [a, b] = pick_line(poly);
        if (!(a == b)) check_cut(poly, a, b, false);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(3, 12);
        long long range = it % 2 ? 5 : 1000;
        vector<Point> pts;
        set<pair<long long, long long>> used;
        while ((int)pts.size() < n){
            Point p = random_point(range);
            if (used.insert({p.x, p.y}).second) pts.push_back(p);
        }
        bool collinear = true;
        for (int i = 2; i < n; i++) collinear &= cross(pts[0], pts[1], pts[i]) == 0;
        if (collinear) continue;

        vector<Point> poly;
        for (int i : simple_polygon(pts)) poly.push_back(pts[i]);
        auto [ex, ey] = brute_centroid(poly);
        long long ox = stress::rand_int(-1000000, 1000000), oy = stress::rand_int(-1000000, 1000000);
        vector<PointF> fpoly, shifted;
        for (Point p : poly) fpoly.push_back(to_f(p)), shifted.push_back({(double)(p.x + ox), (double)(p.y + oy)});
        PointF got = polygon_centroid(fpoly), back = polygon_centroid({fpoly.rbegin(), fpoly.rend()}), far = polygon_centroid(shifted);
        assert(close_to(got.x, ex) && close_to(got.y, ey) && close_to(back.x, ex) && close_to(back.y, ey));
        assert(fabsl(far.x - (ex + ox)) < 1e-6 && fabsl(far.y - (ey + oy)) < 1e-6);
    }

    return 0;
}
