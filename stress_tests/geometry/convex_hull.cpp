#include "../common.h"

#define main library_main
#include "../../code_library/geometry/convex_hull.cpp"
#undef main

bool same(const Point& a, const Point& b){
    return a.x == b.x && a.y == b.y;
}

vector<Point> random_points(int n, int range){
    vector<Point> points;
    for (int i = 0; i < n; i++){
        if (stress::rand_int(0, 4) == 0 && i) points.push_back(points[stress::rand_int(0, i - 1)]);              /// duplicates
        else if (stress::rand_int(0, 5) == 0) points.push_back(Point(stress::rand_int(-range, range), 7));    /// collinear row
        else points.push_back(Point(stress::rand_int(-range, range), stress::rand_int(-range, range)));
    }
    return points;
}

int random_range(){
    return stress::rand_int(0, 2) ? 4 : 1000000000;
}

__int128 cross128(const Point& o, const Point& a, const Point& b){
    return (__int128)(a.x - o.x) * (b.y - o.y) - (__int128)(a.y - o.y) * (b.x - o.x);
}

__int128 dist128(const Point& a, const Point& b){
    return (__int128)(a.x - b.x) * (a.x - b.x) + (__int128)(a.y - b.y) * (a.y - b.y);
}

/// Gift wrapping in 128-bit arithmetic, independent of the monotone chain and safe for Minkowski sums at 2e9
vector<Point> jarvis_hull(vector<Point> points){
    sort(points.begin(), points.end());
    points.erase(unique(points.begin(), points.end(), same), points.end());
    if (points.size() <= 1) return points;
    vector<Point> res;
    Point cur = points[0];

    do {
        res.push_back(cur);
        Point next = same(points[0], cur) ? points[1] : points[0];
        for (auto& p : points){
            __int128 c = cross128(cur, next, p);
            if (c < 0 || (c == 0 && dist128(cur, p) > dist128(cur, next))) next = p;
        }
        cur = next;
    } while (!same(cur, points[0]));
    return res;
}

bool brute_contains(const vector<Point>& hull, const Point& q){
    int h = hull.size();
    if (h == 0) return false;
    if (h == 1) return same(hull[0], q);
    if (h == 2){
        return cross(hull[0], hull[1], q) == 0 && min(hull[0].x, hull[1].x) <= q.x && q.x <= max(hull[0].x, hull[1].x)
            && min(hull[0].y, hull[1].y) <= q.y && q.y <= max(hull[0].y, hull[1].y);
    }

    for (int i = 0; i < h; i++){
        if (cross(hull[i], hull[(i + 1) % h], q) < 0) return false;
    }
    return true;
}

/// The line_hull encoding computed from the side of every vertex
array<int, 2> brute_line_hull(const Point& a, const Point& b, const vector<Point>& hull){
    int h = hull.size(), pos = 0, neg = 0;
    vector<int> s(h);
    for (int i = 0; i < h; i++){
        s[i] = sgn(cross(a, b, hull[i]));
        pos += s[i] > 0, neg += s[i] < 0;
    }
    int zero = h - pos - neg;

    if (!zero && (!pos || !neg)) return {-1, -1};
    if (!pos || !neg){
        assert(zero <= 2);
        for (int i = 0; i < h; i++){
            if (zero == 1 && !s[i]) return {i, -1};
            if (zero == 2 && !s[i] && !s[(i + 1) % h]) return {i, i};
        }
        assert(false);
    }

    array<int, 2> res = {-1, -1};
    for (int i = 0; i < h; i++){
        int before = s[(i + h - 1) % h], here = s[i], after = s[(i + 1) % h];
        if ((here > 0 && after < 0) || (here == 0 && before > 0 && after < 0)) res[0] = i;
        if ((here < 0 && after > 0) || (here == 0 && before < 0 && after > 0)) res[1] = i;
    }
    assert(res[0] != -1 && res[1] != -1);
    return res;
}

/// A strict hull is the unique strictly convex counter-clockwise polygon over input points that contains every input point
void check_hull(const vector<Point>& points, const vector<Point>& hull){
    vector<Point> distinct = points;
    sort(distinct.begin(), distinct.end());
    distinct.erase(unique(distinct.begin(), distinct.end(), same), distinct.end());
    int h = hull.size();

    if (distinct.size() <= 1){
        assert(h == (int)distinct.size() && (h == 0 || same(hull[0], distinct[0])));
        return;
    }

    assert(h >= 2 && same(hull[0], distinct[0]));
    for (auto& v : hull) assert(binary_search(distinct.begin(), distinct.end(), v));
    for (int i = 0; i < h; i++) for (int j = i + 1; j < h; j++) assert(!same(hull[i], hull[j]));

    if (h == 2){  /// every point lies on the segment between the two endpoints
        assert(same(hull[1], distinct.back()));
        for (auto& p : points) assert(cross(hull[0], hull[1], p) == 0);
        return;
    }

    for (int i = 0; i < h; i++){
        assert(cross(hull[i], hull[(i + 1) % h], hull[(i + 2) % h]) > 0);
        for (auto& p : points) assert(cross(hull[i], hull[(i + 1) % h], p) >= 0);
    }
}

bool on_segment(const Point& a, const Point& b, const Point& p){
    return cross(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

/// Distinct vertices form a convex polygon iff all lie on the hull boundary and visit it once around in boundary order
bool brute_is_convex(const vector<Point>& polygon){
    int n = polygon.size();
    auto hull = jarvis_hull(polygon);
    int h = hull.size();
    if (h < 3) return false;

    vector<pair<int, int64_t>> pos(n);
    for (int i = 0; i < n; i++){
        pos[i] = {-1, 0};
        for (int e = 0; e < h && pos[i].first < 0; e++){
            if (on_segment(hull[e], hull[(e + 1) % h], polygon[i]) && !same(polygon[i], hull[(e + 1) % h])) pos[i] = {e, dist2(hull[e], polygon[i])};
        }
        if (pos[i].first < 0) return false;
    }

    int up = 0, down = 0;
    for (int i = 0; i < n; i++) up += pos[i] > pos[(i + 1) % n], down += pos[i] < pos[(i + 1) % n];
    return up == 1 || down == 1;
}

/// Boundary lattice points of a small hull in boundary order, then a swap, a moved vertex, a reversal and a rotation at random
void stress_convexity_orders(){
    for (long long it = 0; it < stress::scaled(100000); it++){
        int range = stress::rand_int(1, 4);
        vector<Point> polygon;

        if (stress::rand_int(0, 1)){
            auto hull = get_convex_hull(random_points(stress::rand_int(3, 8), range));
            int h = hull.size();
            if (h < 3) continue;
            for (int i = 0; i < h; i++){
                Point a = hull[i], b = hull[(i + 1) % h];
                vector<Point> edge;
                for (int x = -range; x <= range; x++) for (int y = -range; y <= range; y++){
                    Point p(x, y);
                    if (on_segment(a, b, p) && !same(p, a) && !same(p, b) && stress::rand_int(0, 1)) edge.push_back(p);
                }
                sort(edge.begin(), edge.end(), [&](const Point& p, const Point& q){ return dist2(a, p) < dist2(a, q); });
                polygon.push_back(a);
                polygon.insert(polygon.end(), edge.begin(), edge.end());
            }
            int m = polygon.size();
            if (stress::rand_int(0, 2) == 0) swap(polygon[stress::rand_int(0, m - 1)], polygon[stress::rand_int(0, m - 1)]);
            if (stress::rand_int(0, 3) == 0){
                Point p(stress::rand_int(-range, range), stress::rand_int(-range, range));
                if (find(polygon.begin(), polygon.end(), p) == polygon.end()) polygon[stress::rand_int(0, m - 1)] = p;
            }
        } else {
            for (int i = stress::rand_int(3, 8); i > 0; i--){
                Point p(stress::rand_int(-range, range), stress::rand_int(-range, range));
                if (find(polygon.begin(), polygon.end(), p) == polygon.end()) polygon.push_back(p);
            }
        }

        if (stress::rand_int(0, 1)) reverse(polygon.begin(), polygon.end());
        if (!polygon.empty()) rotate(polygon.begin(), polygon.begin() + stress::rand_int(0, polygon.size() - 1), polygon.end());
        assert(is_convex(polygon) == brute_is_convex(polygon));
    }
}

void stress_hull_and_convexity(){
    /// Folds along an axis-parallel edge, where only one coordinate changes direction
    assert(!is_convex({Point(0, 0), Point(3, 0), Point(1, 0), Point(4, 0), Point(0, 4)}));
    assert(!is_convex({Point(0, 0), Point(4, 0), Point(0, 4), Point(0, 1), Point(0, 3)}));

    /// A horizontal and a vertical fold whose 180 degree turns cancel, each axis still reverses exactly twice
    assert(!is_convex({Point(0, -1), Point(-1, 0), Point(1, 0), Point(0, 0), Point(0, 1)}));

    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 300), range = random_range();
        vector<Point> points = random_points(n, range);

        auto hull = get_convex_hull(points);
        check_hull(points, hull);

        int h = hull.size();
        if (h < 3){
            assert(!is_convex(hull));
            continue;
        }

        vector<Point> rotated(h);
        int shift = stress::rand_int(0, h - 1);
        for (int i = 0; i < h; i++) rotated[i] = hull[(i + shift) % h];
        assert(is_convex(rotated) && is_convex(vector<Point>(rotated.rbegin(), rotated.rend())));

        /// Inserting the midpoint of an edge keeps it convex, when the midpoint is a lattice point
        Point a = hull[0], b = hull[1];
        if ((a.x + b.x) % 2 == 0 && (a.y + b.y) % 2 == 0){
            vector<Point> with_mid = hull;
            with_mid.insert(with_mid.begin() + 1, Point((a.x + b.x) / 2, (a.y + b.y) / 2));
            assert(is_convex(with_mid));
        }

        if (range == 4){  /// an edge doubling back on itself, a -> 2/3 -> 1/3 -> b, turns one way but is not a simple polygon
            vector<Point> folded;
            for (auto& p : hull) folded.push_back(Point(3 * p.x, 3 * p.y));
            Point s = folded[0], t = folded[1], d((t.x - s.x) / 3, (t.y - s.y) / 3);
            folded.insert(folded.begin() + 1, {Point(s.x + 2 * d.x, s.y + 2 * d.y), Point(s.x + d.x, s.y + d.y)});
            assert(!is_convex(folded));
        }

        if (h >= 4){  /// swapping two adjacent vertices makes the boundary cross itself
            vector<Point> crossed = hull;
            int i = stress::rand_int(0, h - 1);
            swap(crossed[i], crossed[(i + 1) % h]);
            assert(!is_convex(crossed));
        }

        if (h >= 5){  /// visiting every k-th vertex for k coprime with h winds around more than once: a star
            for (int k = 2; k < h - 1; k++){
                if (__gcd(k, h) != 1) continue;
                vector<Point> star;
                for (int i = 0; i < h; i++) star.push_back(hull[(long long)i * k % h]);
                assert(!is_convex(star));
            }
        }
    }
}

/// Farthest pair against all O(n^2) input pairs, width against the slab perpendicular to every input pair, O(n^3)
void stress_calipers(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, it % 10 ? 10 : 30);
        vector<Point> points = random_points(n, random_range());
        auto hull = get_convex_hull(points);

        int64_t best = 0;
        for (auto& p : points) for (auto& q : points) best = max(best, dist2(p, q));
        auto far = farthest_pair(hull);
        assert(dist2(far[0], far[1]) == best);
        assert(find(hull.begin(), hull.end(), far[0]) != hull.end() && find(hull.begin(), hull.end(), far[1]) != hull.end());

        long double slab = numeric_limits<long double>::infinity();
        for (auto& p : points){
            for (auto& q : points){
                if (same(p, q)) continue;
                int64_t lo = 0, hi = 0;
                for (auto& r : points) lo = min(lo, cross(p, q, r)), hi = max(hi, cross(p, q, r));
                slab = min(slab, (hi - lo) / sqrtl(dist2(p, q)));
            }
        }
        if (hull.size() <= 2) slab = 0;
        long double w = width(hull);
        assert(abs(w - slab) <= 1e-9L * max((long double)1, slab));
    }
}

/// Against gift wrapping over every pairwise sum of the input points, inputs rotated and padded with collinear vertices
void stress_minkowski(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int range = random_range();
        vector<Point> a = random_points(stress::rand_int(0, it % 10 ? 8 : 40), range), b = random_points(stress::rand_int(0, it % 10 ? 8 : 40), range);
        vector<Point> ha = get_convex_hull(a), hb = get_convex_hull(b);

        for (auto* poly : {&ha, &hb}){
            int h = poly->size();
            if (h >= 3 && stress::rand_int(0, 1)){
                vector<Point> padded;
                for (int i = 0; i < h; i++){
                    Point p = (*poly)[i], q = (*poly)[(i + 1) % h];
                    padded.push_back(p);
                    if ((p.x + q.x) % 2 == 0 && (p.y + q.y) % 2 == 0) padded.push_back(Point((p.x + q.x) / 2, (p.y + q.y) / 2));
                }
                *poly = padded;
            }
            if (!poly->empty()) rotate(poly->begin(), poly->begin() + stress::rand_int(0, poly->size() - 1), poly->end());
        }

        vector<Point> sums;
        for (auto& p : a) for (auto& q : b) sums.push_back(p + q);
        assert(minkowski_sum(ha, hb) == jarvis_hull(sums));
    }
}

/// Extreme vertex and tangents against a scan of every vertex, line_hull against the encoding derived from vertex sides
void stress_queries(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int range = random_range();
        auto hull = get_convex_hull(random_points(stress::rand_int(1, it % 10 ? 12 : 100), range));
        int h = hull.size();
        if (h == 0) continue;

        int64_t span = range == 4 ? 3 : 2000000000;
        Point dir(stress::rand_int(-span, span), stress::rand_int(-span, span));
        int64_t best = LLONG_MIN;
        for (auto& p : hull) best = max(best, dot(p, dir));
        assert(dot(hull[extreme_vertex(hull, dir)], dir) == best);

        int64_t reach = range == 4 ? 8 : range;
        Point q(stress::rand_int(-reach, reach), stress::rand_int(-reach, reach));
        if (!brute_contains(hull, q)){
            auto [i, j] = tangents(hull, q);
            for (auto& p : hull) assert(cross(q, hull[i], p) >= 0 && cross(q, hull[j], p) <= 0);
        }

        if (h < 3) continue;
        Point a(stress::rand_int(-reach, reach), stress::rand_int(-reach, reach)), b(stress::rand_int(-reach, reach), stress::rand_int(-reach, reach));
        if (stress::rand_int(0, 3) == 0){  /// a line through a vertex or along an edge
            int k = stress::rand_int(0, h - 1);
            a = hull[k], b = stress::rand_int(0, 1) ? hull[(k + 1) % h] : b;
            if (stress::rand_int(0, 1)) swap(a, b);
        }
        if (same(a, b)) continue;
        assert(line_hull(a, b, hull) == brute_line_hull(a, b, hull));
    }
}

/// The online hull against a static rebuild after every insertion
void stress_dynamic(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 20 ? 25 : 200), range = random_range();
        vector<Point> points = random_points(n, range), prefix;
        DynamicHull dynamic;
        int64_t reach = range == 4 ? 6 : range;

        for (auto& p : points){
            dynamic.add(p);
            prefix.push_back(p);
            auto hull = get_convex_hull(prefix);
            assert(dynamic.hull() == hull);

            for (int k = 0; k < 4; k++){
                Point q = k == 0 ? prefix[stress::rand_int(0, prefix.size() - 1)] : Point(stress::rand_int(-reach, reach), stress::rand_int(-reach, reach));
                assert(dynamic.contains(q) == brute_contains(hull, q));
            }
        }
    }
}

/// A 30000 vertex hull on a parabola: linear time queries would need ~1e10 steps, the answers are checked locally
void stress_large(){
    const int n = 30000;
    vector<Point> points;
    for (int i = 0; i < n; i++) points.push_back(Point(i, (int64_t)i * i - 450000000));
    shuffle(points.begin(), points.end(), stress::rng());

    DynamicHull dynamic;
    for (auto& p : points) dynamic.add(p);
    auto hull = get_convex_hull(points);
    assert((int)hull.size() == n && dynamic.hull() == hull);
    auto at = [&](int i){ return hull[((i % n) + n) % n]; };

    for (int it = 0; it < 200000; it++){
        Point dir(stress::rand_int(-2000000000, 2000000000), stress::rand_int(-2000000000, 2000000000));
        int e = extreme_vertex(hull, dir);
        assert(dot(at(e), dir) >= dot(at(e - 1), dir) && dot(at(e), dir) >= dot(at(e + 1), dir));

        Point q = stress::rand_int(0, 1) ? Point(stress::rand_int(-1000000000, 1000000000), stress::rand_int(-1000000000, -450000001))
                                         : Point(stress::rand_int(0, 1) ? -stress::rand_int(1, 1000000000) : n - 1 + stress::rand_int(1, 1000000000), stress::rand_int(-1000000000, 1000000000));
        auto [i, j] = tangents(hull, q);
        assert(cross(q, at(i), at(i - 1)) >= 0 && cross(q, at(i), at(i + 1)) >= 0);
        assert(cross(q, at(j), at(j - 1)) <= 0 && cross(q, at(j), at(j + 1)) <= 0);

        Point c(stress::rand_int(0, n - 1), stress::rand_int(-450000000, 450000000));
        assert(dynamic.contains(c) == (c.y >= c.x * c.x - 450000000 && cross(hull.back(), hull[0], c) >= 0));
    }

    for (int it = 0; it < 300; it++){
        Point a(stress::rand_int(-1000, n + 1000), stress::rand_int(-500000000, 500000000)), b(stress::rand_int(-1000, n + 1000), stress::rand_int(-500000000, 500000000));
        if (!same(a, b)) assert(line_hull(a, b, hull) == brute_line_hull(a, b, hull));
    }
}

int main(){
    stress_hull_and_convexity();
    stress_convexity_orders();
    stress_calipers();
    stress_minkowski();
    stress_queries();
    stress_dynamic();
    stress_large();

    return 0;
}
