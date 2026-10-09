// LINK: -lquadmath
#include "../common.h"
#include <quadmath.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

long double dist(const Point& a, const Point& b){
    return hypotl(a.x - b.x, a.y - b.y);
}

bool close(long double a, long double b, long double scale){
    return fabsl(a - b) <= 1e-9L * max(1.0L, scale);
}

Point rand_point(int range){
    return {(long double)stress::rand_int(-range, range), (long double)stress::rand_int(-range, range)};
}

template<typename G>
long double adaptive_simpson(const G& g, long double a, long double b, long double fa, long double fm, long double fb, long double whole, long double eps, int depth){
    long double m = (a + b) / 2, flm = g((a + m) / 2), frm = g((m + b) / 2);
    long double left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    if (depth == 0 || fabsl(left + right - whole) <= 15 * eps) return left + right + (left + right - whole) / 15;
    return adaptive_simpson(g, a, m, fa, flm, fm, left, eps / 2, depth - 1) + adaptive_simpson(g, m, b, fm, frm, fb, right, eps / 2, depth - 1);
}

/// Adaptive Simpson slab by slab between breakpoints, x = a + (b - a)(1 - cos t) / 2 smooths the square roots at slab ends,
/// adaptivity covers a square root singularity sitting just outside a narrow slab
template<typename F>
long double integrate(F f, vector<long double> xs, long double eps){
    sort(xs.begin(), xs.end());
    long double res = 0;
    for (size_t k = 0; k + 1 < xs.size(); k++){
        long double a = xs[k], b = xs[k + 1];
        if (b - a < 1e-12L) continue;

        auto g = [&](long double t){ return f(a + (b - a) * (1 - cosl(t)) / 2) * (b - a) / 2 * sinl(t); };
        const int pieces = 8;
        for (int s = 0; s < pieces; s++){
            long double lo = PI * s / pieces, hi = PI * (s + 1) / pieces, flo = g(lo), fmid = g((lo + hi) / 2), fhi = g(hi);
            res += adaptive_simpson(g, lo, hi, flo, fmid, fhi, (hi - lo) / 6 * (flo + 4 * fmid + fhi), eps / xs.size() / pieces, 40);
        }
    }

    return res;
}

/// x coordinates where two circle boundaries cross, used only as integration breakpoints
vector<long double> crossing_xs(const Circle& a, const Circle& b){
    long double dx = b.center.x - a.center.x, dy = b.center.y - a.center.y, d = hypotl(dx, dy);
    if (a.r == 0 || b.r == 0 || d == 0 || d > a.r + b.r || d < fabsl(a.r - b.r)) return {};

    long double base = atan2l(dy, dx), half = acosl(clamp((a.r * a.r + d * d - b.r * b.r) / (2 * a.r * d), -1.0L, 1.0L));
    return {a.center.x + a.r * cosl(base + half), a.center.x + a.r * cosl(base - half)};
}

/// Length of the vertical line at x inside the union of the circles, by merging their chords
long double union_chord(const vector<Circle>& circles, long double x){
    vector<pair<long double, long double>> chords;
    for (const Circle& c : circles){
        long double dx = x - c.center.x;
        if (fabsl(dx) >= c.r) continue;
        long double h = sqrtl(c.r * c.r - dx * dx);
        chords.push_back({c.center.y - h, c.center.y + h});
    }
    sort(chords.begin(), chords.end());

    long double res = 0, top = -1e30L;
    for (auto [lo, hi] : chords){
        lo = max(lo, top);
        if (hi > lo) res += hi - lo, top = hi;
    }
    return res;
}

long double union_reference(const vector<Circle>& circles){
    vector<long double> xs;
    long double total = 0;
    for (size_t i = 0; i < circles.size(); i++){
        total += PI * circles[i].r * circles[i].r;
        xs.push_back(circles[i].center.x - circles[i].r);
        xs.push_back(circles[i].center.x + circles[i].r);
        for (size_t j = 0; j < i; j++){
            for (long double x : crossing_xs(circles[i], circles[j])) xs.push_back(x);
        }
    }
    return integrate([&](long double x){ return union_chord(circles, x); }, xs, 1e-9L * max(1.0L, total));
}

/// Length of the vertical line at x inside both the polygon (even-odd crossings) and the circle
long double polygon_chord(const vector<Point>& poly, const Circle& c, long double x){
    long double dx = x - c.center.x;
    if (fabsl(dx) >= c.r) return 0;
    long double h = sqrtl(c.r * c.r - dx * dx), lo = c.center.y - h, hi = c.center.y + h;

    vector<long double> ys;
    for (size_t i = 0; i < poly.size(); i++){
        const Point& p = poly[i];
        const Point& q = poly[(i + 1) % poly.size()];
        if ((p.x < x) != (q.x < x)) ys.push_back(p.y + (q.y - p.y) * (x - p.x) / (q.x - p.x));
    }
    sort(ys.begin(), ys.end());

    long double res = 0;
    for (size_t k = 0; k + 1 < ys.size(); k += 2) res += max(0.0L, min(hi, ys[k + 1]) - max(lo, ys[k]));
    return res;
}

long double polygon_reference(const Circle& c, const vector<Point>& poly){
    vector<long double> xs = {c.center.x - c.r, c.center.x + c.r};
    for (size_t i = 0; i < poly.size(); i++){
        const Point& p = poly[i];
        const Point& q = poly[(i + 1) % poly.size()];
        xs.push_back(p.x);

        /// |p + t (q - p) - center|^2 = r^2 as a quadratic in t
        long double ux = q.x - p.x, uy = q.y - p.y, wx = p.x - c.center.x, wy = p.y - c.center.y;
        long double qa = ux * ux + uy * uy, qb = 2 * (ux * wx + uy * wy), qc = wx * wx + wy * wy - c.r * c.r, disc = qb * qb - 4 * qa * qc;
        if (qa == 0 || disc < 0) continue;
        for (long double t : {(-qb - sqrtl(disc)) / (2 * qa), (-qb + sqrtl(disc)) / (2 * qa)}){
            if (t >= 0 && t <= 1) xs.push_back(p.x + t * ux);
        }
    }
    return integrate([&](long double x){ return polygon_chord(poly, c, x); }, xs, 1e-9L * max(1.0L, c.r * c.r));
}

/// Star-shaped polygon around a non-lattice point, so it is simple and usually not convex
vector<Point> random_polygon(int n, int range){
    while (true){
        vector<Point> poly(n);
        for (auto& p : poly) p = rand_point(range);
        long double ox = 0.1234L, oy = 0.0567L;
        for (auto& p : poly) ox += p.x / n, oy += p.y / n;

        vector<pair<long double, int>> order;
        for (int i = 0; i < n; i++) order.push_back({atan2l(poly[i].y - oy, poly[i].x - ox), i});
        sort(order.begin(), order.end());
        bool distinct = true;
        for (int i = 0; i + 1 < n; i++) distinct &= order[i + 1].first - order[i].first > 1e-9L;
        if (!distinct) continue;

        vector<Point> res;
        for (auto [angle, i] : order) res.push_back(poly[i]);
        if (stress::rand_int(0, 1)) reverse(res.begin(), res.end());
        return res;
    }
}

/// Lens area from the segment formula with the kite area taken by Heron's formula, a different route than the library's sines
__float128 lens_area(__float128 x0, __float128 y0, __float128 r0, __float128 x1, __float128 y1, __float128 r1){
    __float128 d = sqrtq((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)), pi = acosq(-1);
    if (d >= r0 + r1) return 0;
    if (d <= fabsq(r0 - r1)) return pi * fminq(r0, r1) * fminq(r0, r1);

    __float128 a0 = acosq((d * d + r0 * r0 - r1 * r1) / (2 * d * r0)), a1 = acosq((d * d + r1 * r1 - r0 * r0) / (2 * d * r1));
    __float128 kite = sqrtq((-d + r0 + r1) * (d + r0 - r1) * (d - r0 + r1) * (d + r0 + r1)) / 2;
    return r0 * r0 * a0 + r1 * r1 * a1 - kite;
}

/// Center of the circle through three non-collinear points by Cramer's rule on the two perpendicular bisectors
Point cramer_center(const Point& a, const Point& b, const Point& c){
    long double a1 = 2 * (b.x - a.x), b1 = 2 * (b.y - a.y), c1 = b.norm2() - a.norm2();
    long double a2 = 2 * (c.x - a.x), b2 = 2 * (c.y - a.y), c2 = c.norm2() - a.norm2(), det = a1 * b2 - a2 * b1;
    return {(c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det};
}

bool encloses(const vector<Point>& points, const Point& c, long double r){
    for (const Point& p : points){
        if (dist(p, c) > r + 1e-9L * max(1.0L, r)) return false;
    }
    return true;
}

/// Smallest enclosing circle among all single points, pairs as diameters and circumcircles of triples, O(n^4)
long double brute_mec_radius(const vector<Point>& p){
    int n = p.size();
    long double best = 1e30L;
    for (int i = 0; i < n; i++){
        if (encloses(p, p[i], 0)) best = 0;
        for (int j = 0; j < i; j++){
            Point mid{(p[i].x + p[j].x) / 2, (p[i].y + p[j].y) / 2};
            if (encloses(p, mid, dist(p[i], p[j]) / 2)) best = min(best, dist(p[i], p[j]) / 2);
            for (int k = 0; k < j; k++){
                if ((p[j] - p[i]).cross(p[k] - p[i]) == 0) continue;
                /// Relative to p[i], so the squared norms of a tight cluster far from the origin do not cancel
                Point c = p[i] + cramer_center({0, 0}, p[j] - p[i], p[k] - p[i]);
                if (encloses(p, c, dist(c, p[i]))) best = min(best, dist(c, p[i]));
            }
        }
    }
    return best;
}

/// An enclosing circle is the minimum one iff its center lies in the hull of the points on it: no angular gap above pi
bool minimal_certificate(const vector<Point>& points, const Circle& c){
    if (c.r == 0) return true;

    vector<long double> angles;
    for (const Point& p : points){
        if (fabsl(dist(p, c.center) - c.r) <= 1e-9L * c.r) angles.push_back(atan2l(p.y - c.center.y, p.x - c.center.x));
    }
    if (angles.empty()) return false;
    sort(angles.begin(), angles.end());

    long double gap = angles[0] + 2 * PI - angles.back();
    for (size_t i = 0; i + 1 < angles.size(); i++) gap = max(gap, angles[i + 1] - angles[i]);
    return gap <= PI + 1e-7L;
}

int covered_by(const vector<Point>& points, const Point& c, long double r){
    int res = 0;
    for (const Point& p : points) res += (p - c).norm2() <= r * r + 1e-9L * max(1.0L, r * r);
    return res;
}

/// Some optimal circle has two points on its boundary or all its points coincide, so these centers suffice, O(n^3)
int brute_cover(const vector<Point>& points, long double r){
    int n = points.size(), best = 0;
    for (int i = 0; i < n; i++){
        best = max(best, covered_by(points, points[i], r));
        for (int j = 0; j < i; j++){
            Point d = points[j] - points[i];
            long double d2 = d.norm2();
            if (d2 == 0 || d2 > 4 * r * r) continue;
            Point mid = (points[i] + points[j]) / 2, h = d.perp() * (sqrtl(max(0.0L, r * r - d2 / 4)) / sqrtl(d2));
            best = max(best, covered_by(points, mid + h, r));
            best = max(best, covered_by(points, mid - h, r));
        }
    }
    return best;
}

const int TRIPLES[][3] = {{1, 0, 1}, {0, 1, 1}, {3, 4, 5}, {4, 3, 5}, {5, 12, 13}, {8, 15, 17}, {7, 24, 25}};

/// A random unit direction with rational coordinates dx / len, dy / len, so tangent configurations stay exact on integers,
/// the checks scale some of them by 0.1, which no binary float holds, so the tolerances decide those tangencies
array<long long, 3> rational_direction(){
    auto t = TRIPLES[stress::rand_int(0, 6)];
    long long sx = stress::rand_int(0, 1) ? 1 : -1, sy = stress::rand_int(0, 1) ? 1 : -1;
    return {t[0] * sx, t[1] * sy, t[2]};
}

vector<Point> ring_points(long long r){
    vector<Point> res;
    for (long long x = -r; x <= r; x++){
        long long y2 = r * r - x * x, y = llround(sqrtl(y2));
        if (y * y != y2) continue;
        res.push_back({(long double)x, (long double)y});
        if (y != 0) res.push_back({(long double)x, (long double)-y});
    }
    return res;
}

int sign128(__int128 v){
    return (v > 0) - (v < 0);
}

void check_line_intersections(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int range = it % 3 ? 20 : 10000;
        long double unit = 1;
        if (it % 8 < 2) range *= 10, unit = 0.1L;
        /// A zero radius has zero tolerance, so it only meets a line exactly and is left to the integer iterations
        Point c = rand_point(range), a, b;
        long long r = stress::rand_int(unit < 1, range);
        if (it % 4 == 0){
            auto [dx, dy, len] = rational_direction();
            long long k = stress::rand_int(unit < 1, range / 25 + 1), m = stress::rand_int(-5, 5);
            r = k * len;
            a = {c.x + dx * k + dy * m, c.y + dy * k - dx * m};
            b = {a.x + dy, a.y - dx};
        }
        else{
            do a = rand_point(range), b = rand_point(range); while (a.x == b.x && a.y == b.y);
        }
        if (it % 2) swap(a, b);

        __int128 ux = (long long)(b.x - a.x), uy = (long long)(b.y - a.y), wx = (long long)(c.x - a.x), wy = (long long)(c.y - a.y);
        __int128 s = ux * wy - uy * wx, len2 = ux * ux + uy * uy;
        int expected = 1 + sign128((__int128)r * r * len2 - s * s);

        c = c * unit, a = a * unit, b = b * unit;
        long double sr = r * unit;
        vector<Point> res = circle_line_intersection({c, sr}, a, b);
        assert((int)res.size() == expected);
        for (const Point& p : res){
            assert(close(dist(p, c), sr, sr));
            assert(close((b - a).cross(p - a) / (b - a).norm(), 0, range * unit));
        }
        if (expected == 2) assert((res[1] - res[0]).dot(b - a) > 0);
    }
}

void check_circle_intersections(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int range = it % 3 ? 20 : 10000;
        long double unit = 1;
        if (it % 8 < 2) range *= 10, unit = 0.1L;
        Point ca = rand_point(range), cb = rand_point(range);
        long long ra = stress::rand_int(0, range), rb = stress::rand_int(0, range);
        if (it % 4 == 0){
            auto [dx, dy, len] = rational_direction();
            long long k = stress::rand_int(1, range / 25 + 1);
            cb = {ca.x + dx * k, ca.y + dy * k};
            ra = stress::rand_int(0, len * k);
            rb = stress::rand_int(0, 1) ? len * k - ra : len * k + ra;
        }
        if (it % 50 == 0) cb = ca;

        long long dx = cb.x - ca.x, dy = cb.y - ca.y, d2 = dx * dx + dy * dy;
        int expected = 0;
        if (d2 != 0 && d2 <= (ra + rb) * (ra + rb) && d2 >= (ra - rb) * (ra - rb)){
            expected = d2 == (ra + rb) * (ra + rb) || d2 == (ra - rb) * (ra - rb) ? 1 : 2;
        }

        ca = ca * unit, cb = cb * unit;
        long double sa = ra * unit, sb = rb * unit;
        vector<Point> res = circle_circle_intersection({ca, sa}, {cb, sb});
        assert((int)res.size() == expected);
        for (const Point& p : res){
            assert(close(dist(p, ca), sa, max(sa, sb)));
            assert(close(dist(p, cb), sb, max(sa, sb)));
        }
        if (expected == 2) assert(dist(res[0], res[1]) > 0);
    }
}

void check_point_tangents(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int range = it % 3 ? 20 : 10000;
        long double unit = 1;
        if (it % 8 < 2) range *= 10, unit = 0.1L;
        Point c = rand_point(range), p = rand_point(range);
        long long r = stress::rand_int(1, range);
        if (it % 4 == 0){
            auto [dx, dy, len] = rational_direction();
            long long k = stress::rand_int(1, range / 25 + 1);
            r = len * k, p = {c.x + dx * k, c.y + dy * k};
        }

        long long dx = p.x - c.x, dy = p.y - c.y, d2 = dx * dx + dy * dy;
        int expected = 1 + (d2 > r * r) - (d2 < r * r);

        c = c * unit, p = p * unit;
        long double sr = r * unit;
        vector<Point> res = tangents_from_point({c, sr}, p);
        assert((int)res.size() == expected);
        for (const Point& t : res){
            assert(close(dist(t, c), sr, sr));
            assert(close((t - c).dot(t - p), 0, sr * sqrtl(d2) * unit));
        }
        if (expected == 1) assert(dist(res[0], p) == 0);
        if (expected == 2) assert(dist(res[0], res[1]) > 1e-9L * sr);
    }
}

void check_common_tangents(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int range = it % 3 ? 20 : 10000;
        long double unit = 1;
        if (it % 8 < 2) range *= 10, unit = 0.1L;
        Point ca = rand_point(range), cb = rand_point(range);
        long long ra = stress::rand_int(1, range), rb = stress::rand_int(1, range);
        if (it % 4 == 0){
            auto [dx, dy, len] = rational_direction();
            long long k = stress::rand_int(1, range / 25 + 1);
            cb = {ca.x + dx * k, ca.y + dy * k};
            ra = stress::rand_int(1, len * k);
            rb = stress::rand_int(0, 1) && len * k > ra ? len * k - ra : len * k + ra;
        }
        if (it % 50 == 0) cb = ca;
        if (it % 100 == 0) rb = ra;

        long long dx = cb.x - ca.x, dy = cb.y - ca.y, d2 = dx * dx + dy * dy;
        auto lines = [&](long long gap){ return d2 == 0 ? 0 : 1 + (d2 > gap * gap) - (d2 < gap * gap); };
        int outer = lines(ra - rb), inner = lines(ra + rb);

        Circle a{ca * unit, ra * unit}, b{cb * unit, rb * unit};
        vector<pair<Point, Point>> res = common_tangents(a, b);
        assert((int)res.size() == outer + inner);
        for (int i = 0; i < (int)res.size(); i++){
            auto [ta, tb] = res[i];
            long double scale = (max(ra, rb) + sqrtl(d2)) * unit;
            assert(close(dist(ta, a.center), a.r, scale));
            assert(close(dist(tb, b.center), b.r, scale));
            for (int j = 0; j < i; j++) assert(dist(res[j].first, ta) > 1e-9L * a.r);

            if (dist(ta, tb) <= 1e-9L * scale) continue;
            Point u = (tb - ta) / dist(ta, tb);
            assert(close(u.dot(ta - a.center), 0, scale));
            assert(close(u.dot(tb - b.center), 0, scale));

            bool same_side = (u.cross(a.center - ta) > 0) == (u.cross(b.center - ta) > 0);
            assert(same_side == (i < outer));
        }
    }
}

void check_circumcircle(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int range = it % 3 ? 20 : 10000;
        Point a = rand_point(range), b = rand_point(range), c = rand_point(range);
        if ((b - a).cross(c - a) == 0) continue;

        Circle res = circumcircle(a, b, c);
        Point expected = cramer_center(a, b, c);
        long double scale = dist(expected, a);
        assert(close(res.center.x, expected.x, scale) && close(res.center.y, expected.y, scale));
        assert(close(res.r, dist(expected, a), scale));
        assert(close(dist(res.center, b), res.r, scale) && close(dist(res.center, c), res.r, scale));
    }
}

void check_minimum_enclosing_circle(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int range = it % 3 == 0 ? 2 : (it % 3 == 1 ? 10 : 10000), n = it < 50 ? 1 + it % 3 : stress::rand_int(1, 9);
        vector<Point> points(n);
        for (auto& p : points) p = rand_point(range);
        if (it % 7 == 0) for (auto& p : points) p.y = 2 * p.x + 1;

        Circle res = minimum_enclosing_circle(points);
        assert(encloses(points, res.center, res.r));
        assert(close(res.r, brute_mec_radius(points), range));
    }

    /// Clusters about 1e-6 wide near coordinates up to 1e4, where rounding the center exceeds the relative tolerance of a tiny radius,
    /// sometimes with farther points so the cluster decides only an intermediate circle
    for (long long it = 0; it < stress::scaled(20000); it++){
        Point base = rand_point(10000);
        int n = stress::rand_int(2, 8), far = it % 2 ? stress::rand_int(1, 3) : 0;
        vector<Point> points;
        for (int i = 0; i < n; i++){
            if (i && stress::rand_int(0, 2) == 0) points.push_back(points[stress::rand_int(0, i - 1)]);
            else points.push_back({base.x + stress::rand_int(0, 20) / 1e7L, base.y + stress::rand_int(0, 20) / 1e7L});
        }
        for (int i = 0; i < far; i++) points.push_back(base + rand_point(100));

        Circle res = minimum_enclosing_circle(points);
        assert(encloses(points, res.center, res.r));
        assert(close(res.r, brute_mec_radius(points), 10000));
    }

    /// Large inputs for the expected O(n) bound: random squares, and lattice points on one circle where every point is on the boundary
    for (int round = 0; round < 6; round++){
        int n = stress::rand_int(50000, 100000);
        vector<Point> points;
        if (round % 2 == 0){
            for (int i = 0; i < n; i++) points.push_back(rand_point(10000));
        }
        else{
            points = ring_points(5525);
            int base = points.size();
            for (int i = 0; i < n; i++) points.push_back(points[stress::rand_int(0, base - 1)]);
        }

        Circle res = minimum_enclosing_circle(points);
        assert(encloses(points, res.center, res.r));
        assert(minimal_certificate(points, res));
        if (round % 2) assert(close(res.r, 5525, 5525) && close(res.center.x, 0, 5525) && close(res.center.y, 0, 5525));
    }
}

void check_polygon_area(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int range = it % 2 ? 6 : 10000, n = stress::rand_int(3, 9);
        vector<Point> poly = random_polygon(n, range);
        Circle c{rand_point(range), (long double)stress::rand_int(0, range)};
        if (it % 5 == 0) c.r += stress::rand_int(0, 999) / 1000.0L;
        if (it % 7 == 0) c.center = poly[stress::rand_int(0, n - 1)];

        long double expected = polygon_reference(c, poly), scale = max(1.0L, c.r * c.r);
        assert(close(circle_polygon_area(c, poly), expected, scale * 100));

        int at = stress::rand_int(0, n - 1);
        poly.insert(poly.begin() + at, poly[at]);
        assert(close(circle_polygon_area(c, poly), expected, scale * 100));
    }
}

void check_circle_circle_area(){
    /// acos(-1) in double is about 1e-4 off for a contained circle of radius 1e6
    assert(fabsl(PI - (long double)acosq(-1)) <= 1e-18L);

    Circle huge{{0, 0}, 2e6}, inner{{1, 1}, 1e6};
    assert(fabsl(circle_circle_area(huge, inner) - (long double)lens_area(0, 0, 2e6, 1, 1, 1e6)) <= 1e-6L);

    for (long long it = 0; it < stress::scaled(200000); it++){
        int range = it % 3 ? 10 : 10000;
        long double x0 = stress::rand_int(-range, range), y0 = stress::rand_int(-range, range), r0 = stress::rand_int(0, range);
        long double x1 = stress::rand_int(-range, range), y1 = stress::rand_int(-range, range), r1 = stress::rand_int(0, range);
        if (it % 4 == 0) x0 += stress::rand_int(0, 999) / 1000.0L, r1 += stress::rand_int(0, 999) / 1000.0L;

        Circle a{{x0, y0}, r0}, b{{x1, y1}, r1};
        long double res = circle_circle_area(a, b), swapped = circle_circle_area(b, a);
        long double expected = lens_area(x0, y0, r0, x1, y1, r1);
        long double scale = max(1.0L, max(r0, r1) * max(r0, r1));

        assert(!std::isnan(res));
        assert(fabsl(res - expected) <= 1e-9L * scale);
        assert(fabsl(res - swapped) <= 1e-9L * scale);
        assert(res >= -1e-9L * scale && res <= PI * min(r0, r1) * min(r0, r1) + 1e-9L * scale);

        if (it % 500 == 0) assert(fabsl(res - (PI * r0 * r0 + PI * r1 * r1 - union_reference({a, b}))) <= 1e-7L * scale);
    }

    /// Tangent and nearly tangent pairs, where the acos arguments sit right at +-1
    for (long long it = 0; it < stress::scaled(20000); it++){
        long double r0 = stress::rand_int(1, 1000), r1 = stress::rand_int(1, 1000), t = stress::rand_int(-1000, 1000) * 1e-12L;
        long double d = stress::rand_int(0, 1) ? r0 + r1 + t : fabsl(r0 - r1) + t;
        long double angle = stress::rand_int(0, 359) * PI / 180;
        Circle a{{0, 0}, r0}, b{{d * cosl(angle), d * sinl(angle)}, r1};

        long double res = circle_circle_area(a, b), scale = max(r0, r1) * max(r0, r1);
        assert(!std::isnan(res));
        assert(fabsl(res - (long double)lens_area(0, 0, r0, b.center.x, b.center.y, r1)) <= 1e-6L * scale);
    }
}

void check_union_area(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int range = it % 3 ? 5 : 10000, n = it < 20 ? it % 3 : stress::rand_int(1, 7);
        vector<Circle> circles;
        for (int i = 0; i < n; i++){
            if (i > 0 && stress::rand_int(0, 4) == 0){
                circles.push_back(circles[stress::rand_int(0, i - 1)]);
                continue;
            }
            Circle c{rand_point(range), (long double)stress::rand_int(0, range)};
            if (it % 4 == 0) c.center.x += stress::rand_int(0, 999) / 1000.0L, c.r += stress::rand_int(0, 999) / 1000.0L;
            circles.push_back(c);
        }

        long double total = 0;
        for (const Circle& c : circles) total += PI * c.r * c.r;
        long double res = circle_union_area(circles), expected = union_reference(circles);
        assert(fabsl(res - expected) <= 1e-7L * max(1.0L, total));
        if (n == 2) assert(close(res, total - circle_circle_area(circles[0], circles[1]), total));
    }

    for (int round = 0; round < 2; round++){
        vector<Circle> circles;
        for (int i = 0; i < 30; i++) circles.push_back({rand_point(100), (long double)stress::rand_int(1, 30)});

        long double total = 0;
        for (const Circle& c : circles) total += PI * c.r * c.r;
        assert(fabsl(circle_union_area(circles) - union_reference(circles)) <= 1e-7L * total);
    }

    /// n = 2000 guards the O(n^2 log n) bound, the reference is too slow there so only bounds are checked
    vector<Circle> circles;
    long double total = 0, largest = 0;
    for (int i = 0; i < 2000; i++){
        circles.push_back({rand_point(3000), (long double)stress::rand_int(1, 1000)});
        total += PI * circles[i].r * circles[i].r;
        largest = max(largest, PI * circles[i].r * circles[i].r);
    }

    long double res = circle_union_area(circles);
    assert(res >= largest * (1 - 1e-9L) && res <= total * (1 + 1e-9L) && res <= 8000.0L * 8000.0L);
}

void check_max_cover(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int range = it % 4 == 3 ? 10000 : (it % 2 ? 4 : 8), n = it < 20 ? it % 4 : stress::rand_int(1, 14);
        vector<Point> points(n);
        for (auto& p : points) p = rand_point(range);
        long double r = it % 3 == 0 ? stress::rand_int(0, 2 * range) / 2.0L : stress::rand_int(0, 500 * range) / 1000.0L;

        auto [count, center] = max_circle_cover(points, r);
        assert(count == brute_cover(points, r));
        if (n > 0) assert(covered_by(points, center, r) == count);
    }

    for (int round = 0; round < 3; round++){
        vector<Point> points(150);
        for (auto& p : points) p = rand_point(30);
        long double r = stress::rand_int(1, 10);

        auto [count, center] = max_circle_cover(points, r);
        assert(count == brute_cover(points, r) && covered_by(points, center, r) == count);
    }

    /// Lattice points on one radius-r circle tie exactly on the boundary of the optimal circle
    const long long RADII[] = {5, 13, 25, 65, 85};
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long r = RADII[stress::rand_int(0, 4)];
        Point offset = rand_point(10000 - 2 * r);
        vector<Point> points;
        for (const Point& p : ring_points(r)){
            if (stress::rand_int(0, 2)) points.push_back(p + offset);
        }

        int on_ring = points.size(), extra = it % 2 ? stress::rand_int(1, 4) : 0;
        for (int i = 0; i < extra; i++) points.push_back(offset + rand_point(2 * r));
        shuffle(points.begin(), points.end(), stress::rng());

        auto [count, center] = max_circle_cover(points, r);
        assert(count == brute_cover(points, r) && covered_by(points, center, r) == count);
        if (extra == 0) assert(count == on_ring);
    }

    /// n = 2000 guards the O(n^2 log n) bound, the brute force is too slow there so only bounds are checked
    vector<Point> points(2000);
    for (auto& p : points) p = rand_point(300);
    long double r = stress::rand_int(30, 100);

    int single = 0;
    for (const Point& p : points) single = max(single, covered_by(points, p, r));
    auto [count, center] = max_circle_cover(points, r);
    assert(count >= single && count <= 2000 && covered_by(points, center, r) == count);
}

/// Intersections, tangents and the circumcircle against exact integer classification and the defining properties,
/// the enclosing circle against all O(n^3) candidates, areas against numerical integration, cover against candidate centers
int main(){
    check_line_intersections();
    check_circle_intersections();
    check_point_tangents();
    check_common_tangents();
    check_circumcircle();
    check_minimum_enclosing_circle();
    check_polygon_area();
    check_circle_circle_area();
    check_union_area();
    check_max_cover();

    return 0;
}
