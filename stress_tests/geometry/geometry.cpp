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

    return 0;
}
