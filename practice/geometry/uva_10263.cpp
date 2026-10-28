/***
 *
 * 2D Geometry
 * Exact integer predicates on lattice points, polygon queries, and floating point line and polygon helpers
 *
 * Complexity, n = number of polygon vertices or points:
 *   - O(n log n) for simple_polygon
 *   - O(log n) for point_in_convex_polygon
 *   - O(n) for area2, point_in_polygon, cut_polygon, polygon_centroid
 *   - O(1) for every other function
 *
 * Integer part, Point of long long, |x|, |y| <= 1e9 so every cross / dot product fits in long long:
 *   Point has +, -, == and < (x, then y), so points sort, deduplicate and key a map
 *   cross(o, a, b): > 0 if o -> a -> b turns left (counter-clockwise), < 0 right, 0 collinear
 *   on_segment(p, a, b), segments_intersect(a, b, c, d): endpoints inclusive, degenerate segments allowed
 *   dist_point_segment(p, a, b): Euclidean distance as a double
 *   area2(poly): twice the signed area, positive for counter-clockwise order
 *   point_in_polygon(poly, p): any simple polygon, returns 1 inside, 0 on the boundary, -1 outside
 *   point_in_convex_polygon(poly, p): convex polygon in counter-clockwise order with n >= 3, same return values
 *       collinear vertices allowed, as long as not all vertices are collinear
 *   simple_polygon(points): an order of the indices forming a simple polygon
 *       points must be distinct and not all collinear
 *   angle_less(a, b): polar angle order of vectors, counter-clockwise from the positive x axis over [0, 2 pi)
 *       vectors with the same direction compare equal, (0, 0) comes first
 *       components up to 2e9 are allowed, so sort around a center o by sorting the vectors p - o
 *
 * Floating point part, PointF of doubles:
 *   Tolerance is relative to the size of the configuration: o, a, b count as collinear when the height of the
 *   triangle is below EPS = 1e-9 times its longest side, two lines are parallel when the sine of their angle is below EPS
 *   With integer coordinates |x|, |y| <= 1e4 every decision below is exact, only the returned points are rounded
 *   orientation(o, a, b): sign of cross(o, a, b) under that tolerance
 *   line_intersection(a, b, c, d): lines through a, b and c, d with a != b, c != d
 *       returns {1, point} for a unique point, {0, (0, 0)} for parallel lines, {-1, (0, 0)} for the same line
 *   segment_intersection(a, b, c, d): no point, the single common point, or the two ends of the overlap of collinear
 *       segments, degenerate segments allowed
 *   dist_point_line(p, a, b): signed distance, positive when p is left of a -> b, a != b
 *   project(p, a, b), reflect(p, a, b): projection of p onto, reflection of p across the line through a, b, a != b
 *   cut_polygon(poly, a, b): the part of poly left of a -> b, the half-plane cross(a, b, p) >= 0
 *       a != b, decided by the exact sign of cross(a, b, p) without the EPS tolerance
 *       empty when no vertex lies strictly left of the line
 *       convex poly: a convex polygon with the same orientation, no repeated vertex for exact input
 *       (a vertex within rounding of the line may appear twice)
 *       simple poly: the kept pieces joined by zero-width bridges along the line, the signed area is exact
 *   area2(poly): twice the signed area, positive for counter-clockwise order, 0 for fewer than 3 vertices
 *   polygon_centroid(poly): center of mass of a simple polygon with nonzero area, either orientation
 *   rotate(center, p, degrees): rotates p counter-clockwise around center
 *   clockwise_angle(a, b): clockwise angle in degrees in [0, 360) from vector a to vector b
 *   great_circle_distance(lat1, lon1, lat2, lon2, radius): haversine formula, angles in degrees
 *
 * Defines its own Point, cross and area2: pasted next to another geometry file, merge the Point blocks first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;

    Point operator+(const Point& p) const{
        return {x + p.x, y + p.y};
    }

    Point operator-(const Point& p) const{
        return {x - p.x, y - p.y};
    }

    bool operator==(const Point& p) const{
        return x == p.x && y == p.y;
    }

    bool operator<(const Point& p) const{
        return x != p.x ? x < p.x : y < p.y;
    }
};

long long cross(const Point& o, const Point& a, const Point& b){
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

long long dot(const Point& o, const Point& a, const Point& b){
    return (a.x - o.x) * (b.x - o.x) + (a.y - o.y) * (b.y - o.y);
}

bool on_segment(const Point& p, const Point& a, const Point& b){
    return cross(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

bool segments_intersect(const Point& a, const Point& b, const Point& c, const Point& d){
    auto sign = [](long long v){ return (v > 0) - (v < 0); };
    int d1 = sign(cross(a, b, c)), d2 = sign(cross(a, b, d)), d3 = sign(cross(c, d, a)), d4 = sign(cross(c, d, b));
    if (d1 * d2 < 0 && d3 * d4 < 0) return true;
    return on_segment(c, a, b) || on_segment(d, a, b) || on_segment(a, c, d) || on_segment(b, c, d);
}

double dist_point_segment(const Point& p, const Point& a, const Point& b){
    long long len = dot(a, b, b), t = dot(a, b, p);
    if (len == 0 || t <= 0) return hypot((double)(p.x - a.x), (double)(p.y - a.y));
    if (t >= len) return hypot((double)(p.x - b.x), (double)(p.y - b.y));
    return fabs((double)cross(a, b, p)) / sqrt((double)len);
}

long long area2(const vector<Point>& poly){
    __int128 res = 0;  /// partial sums can exceed long long even when the final area fits
    for (int i = 1; i + 1 < (int)poly.size(); i++) res += cross(poly[0], poly[i], poly[i + 1]);
    return (long long)res;
}

int point_in_polygon(const vector<Point>& poly, const Point& p){
    bool inside = false;
    for (int i = 0, n = poly.size(); i < n; i++){
        const Point& a = poly[i];
        const Point& b = poly[(i + 1) % n];
        if (on_segment(p, a, b)) return 0;
        if ((a.y > p.y) != (b.y > p.y)){
            long long c = cross(a, b, p);
            if ((b.y > a.y) ? c > 0 : c < 0) inside = !inside;
        }
    }

    return inside ? 1 : -1;
}

int point_in_convex_polygon(const vector<Point>& poly, const Point& p){
    int n = poly.size();
    assert(n >= 3);
    const Point& o = poly[0];
    long long first = cross(o, poly[1], p), last = cross(o, poly[n - 1], p);
    if (first < 0 || last > 0) return -1;
    if (first == 0 || last == 0){
        /// The first and last edges may run on through collinear vertices, so test against the farthest one on each ray
        auto on_ray = [&](const Point& dir, const Point& v){ return cross(o, dir, v) == 0 && dot(o, dir, v) > 0; };
        int a = 1, hi = n - 1;
        while (a < hi){
            int mid = (a + hi + 1) / 2;
            if (on_ray(poly[1], poly[mid])) a = mid;
            else hi = mid - 1;
        }

        int lo = 1, b = n - 1;
        while (lo < b){
            int mid = (lo + b) / 2;
            if (on_ray(poly[n - 1], poly[mid])) b = mid;
            else lo = mid + 1;
        }

        return on_segment(p, o, poly[a]) || on_segment(p, o, poly[b]) ? 0 : -1;
    }

    int lo = 1, hi = n - 1;  /// p lies strictly inside the angle, find the wedge o, poly[lo], poly[lo + 1]
    while (hi - lo > 1){
        int mid = (lo + hi) / 2;
        if (cross(o, poly[mid], p) > 0) lo = mid;
        else hi = mid;
    }
    long long c = cross(poly[lo], poly[hi], p);
    return c > 0 ? 1 : (c == 0 ? 0 : -1);
}

vector<int> simple_polygon(const vector<Point>& points){
    int n = points.size();
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    int pivot = *min_element(order.begin(), order.end(), [&](int i, int j){
        return make_pair(points[i].y, points[i].x) < make_pair(points[j].y, points[j].x);
    });
    swap(order[0], order[pivot]);

    const Point& o = points[order[0]];
    sort(order.begin() + 1, order.end(), [&](int i, int j){
        long long c = cross(o, points[i], points[j]);
        if (c != 0) return c > 0;
        return dot(o, points[i], points[i]) < dot(o, points[j], points[j]);
    });

    /// Points on the last ray from the pivot are walked back towards it, otherwise the polygon folds over itself
    int k = n - 1;
    while (k > 1 && cross(o, points[order[k]], points[order[k - 1]]) == 0) k--;
    reverse(order.begin() + k, order.end());
    return order;
}

bool angle_less(const Point& a, const Point& b){
    auto half = [](const Point& p){ return p.x == 0 && p.y == 0 ? 0 : 1 + (p.y < 0 || (p.y == 0 && p.x < 0)); };
    int ha = half(a), hb = half(b);
    if (ha != hb) return ha < hb;
    return cross({0, 0}, a, b) > 0;
}

const double EPS = 1e-9;

struct PointF{
    double x, y;

    PointF operator+(const PointF& p) const{
        return {x + p.x, y + p.y};
    }

    PointF operator-(const PointF& p) const{
        return {x - p.x, y - p.y};
    }

    PointF operator*(double k) const{
        return {x * k, y * k};
    }

    PointF operator/(double k) const{
        return {x / k, y / k};
    }
};

double cross(const PointF& o, const PointF& a, const PointF& b){
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

double dot(const PointF& o, const PointF& a, const PointF& b){
    return (a.x - o.x) * (b.x - o.x) + (a.y - o.y) * (b.y - o.y);
}

int orientation(const PointF& o, const PointF& a, const PointF& b){
    double c = cross(o, a, b), tol = EPS * max({dot(o, a, a), dot(o, b, b), dot(a, b, b)});
    return c > tol ? 1 : (c < -tol ? -1 : 0);
}

pair<int, PointF> line_intersection(const PointF& a, const PointF& b, const PointF& c, const PointF& d){
    double den = (b.x - a.x) * (d.y - c.y) - (b.y - a.y) * (d.x - c.x);
    if (fabs(den) <= EPS * sqrt(dot(a, b, b) * dot(c, d, d))) return {orientation(a, b, c) == 0 ? -1 : 0, {0, 0}};
    return {1, a + (b - a) * (cross(a, c, d) / den)};
}

vector<PointF> segment_intersection(const PointF& a, const PointF& b, const PointF& c, const PointF& d){
    int oa = orientation(c, d, a), ob = orientation(c, d, b), oc = orientation(a, b, c), od = orientation(a, b, d);
    if (oa * ob < 0 && oc * od < 0){
        double ca = cross(c, d, a), cb = cross(c, d, b);
        return {(a * cb - b * ca) / (cb - ca)};
    }

    auto touches = [](const PointF& p, const PointF& s, const PointF& e){
        return orientation(s, e, p) == 0 && dot(p, s, e) <= EPS * dot(s, e, e);
    };
    vector<PointF> found;
    if (touches(a, c, d)) found.push_back(a);
    if (touches(b, c, d)) found.push_back(b);
    if (touches(c, a, b)) found.push_back(c);
    if (touches(d, a, b)) found.push_back(d);
    if (found.empty()) return found;

    /// The overlap ends are the farthest pair, lexicographic order is unreliable when the shared line is near vertical
    PointF u = found[0], v = found[0];
    for (const PointF& p : found){
        for (const PointF& q : found){
            if (dot(p, q, q) > dot(u, v, v)) u = p, v = q;
        }
    }
    if (dot(u, v, v) <= EPS * EPS * max(dot(a, b, b), dot(c, d, d))) return {u};
    return {u, v};
}

double dist_point_line(const PointF& p, const PointF& a, const PointF& b){
    return cross(a, b, p) / hypot(b.x - a.x, b.y - a.y);
}

PointF project(const PointF& p, const PointF& a, const PointF& b){
    return a + (b - a) * (dot(a, b, p) / dot(a, b, b));
}

PointF reflect(const PointF& p, const PointF& a, const PointF& b){
    return project(p, a, b) * 2 - p;
}

vector<PointF> cut_polygon(const vector<PointF>& poly, const PointF& a, const PointF& b){
    vector<PointF> res;
    bool any_left = false;
    for (int i = 0, n = poly.size(); i < n; i++){
        const PointF& cur = poly[i];
        const PointF& prev = poly[(i + n - 1) % n];
        double sc = cross(a, b, cur), sp = cross(a, b, prev);
        /// Only strict sign changes are interpolated: a vertex on the line is kept once as itself, never again as a crossing
        if ((sc > 0 && sp < 0) || (sc < 0 && sp > 0)) res.push_back(cur + (prev - cur) * (sc / (sc - sp)));
        if (sc >= 0) res.push_back(cur);
        any_left |= sc > 0;
    }

    if (!any_left) res.clear();  /// only vertices on the line, a zero-area remainder
    return res;
}

double area2(const vector<PointF>& poly){
    double res = 0;  /// fanned from poly[0] like polygon_centroid, so far away polygons keep their precision
    for (int i = 1; i + 1 < (int)poly.size(); i++) res += cross(poly[0], poly[i], poly[i + 1]);
    return res;
}

PointF polygon_centroid(const vector<PointF>& poly){
    /// Fan from poly[0] instead of the origin, so far away polygons do not lose precision to cancellation
    const PointF& o = poly[0];
    PointF sum{0, 0};
    double area = 0;
    for (int i = 1; i + 1 < (int)poly.size(); i++){
        double c = cross(o, poly[i], poly[i + 1]);
        sum = sum + (poly[i] + poly[i + 1] - o * 2) * c;
        area += c;
    }

    return o + sum / (3 * area);
}

PointF rotate(const PointF& center, const PointF& p, double degrees){
    double theta = degrees * acos(-1.0) / 180.0, s = sin(theta), c = cos(theta);
    double x = p.x - center.x, y = p.y - center.y;
    return {center.x + x * c - y * s, center.y + x * s + y * c};
}

double clockwise_angle(const PointF& a, const PointF& b){
    double theta = atan2(b.x * a.y - b.y * a.x, a.x * b.x + a.y * b.y) * 180.0 / acos(-1.0);
    if (theta < 0) theta += 360.0;
    return theta >= 360.0 ? 0.0 : theta;  /// a tiny negative angle rounds up to exactly 360
}

long double great_circle_distance(long double lat1, long double lon1, long double lat2, long double lon2, long double radius){
    const long double to_rad = acosl(-1.0L) / 180.0L;
    lat1 *= to_rad, lat2 *= to_rad, lon1 *= to_rad, lon2 *= to_rad;
    long double h = powl(sinl((lat2 - lat1) / 2), 2) + cosl(lat1) * cosl(lat2) * powl(sinl((lon2 - lon1) / 2), 2);
    return 2 * radius * asinl(sqrtl(min(1.0L, h)));
}

int main(){
    PointF m;

    while (scanf("%lf %lf", &m.x, &m.y) == 2){
        int n;
        scanf("%d", &n);
        vector<PointF> line(n + 1);
        for (auto& p: line) scanf("%lf %lf", &p.x, &p.y);

        PointF best = line[0];
        double best_dist = dot(m, best, best);
        for (int i = 0; i < n; i++){
            PointF a = line[i], b = line[i + 1], c = a;
            if (dot(a, b, b) > 0 && dot(a, b, m) > 0) c = dot(b, a, m) > 0 ? project(m, a, b) : b;

            double d = dot(m, c, c);
            if (d < best_dist - 1e-9 * max(1.0, best_dist)) best = c, best_dist = d;
        }
        printf("%.4f\n%.4f\n", best.x, best.y);
    }

    return 0;
}
