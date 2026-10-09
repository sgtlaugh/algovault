/***
 *
 * 2D Geometry
 * Exact integer predicates on lattice points, polygon queries, and a few floating point helpers
 *
 * Integer part, |x|, |y| <= 1e9 so every cross / dot product fits in long long:
 *   cross(o, a, b): > 0 if o -> a -> b turns left (counter-clockwise), < 0 right, 0 collinear
 *   on_segment(p, a, b), segments_intersect(a, b, c, d): endpoints inclusive, degenerate segments allowed
 *   dist_point_segment(p, a, b): Euclidean distance as a double
 *   area2(poly): twice the signed area, positive for counter-clockwise order
 *   point_in_polygon(poly, p): any simple polygon, O(n), returns 1 inside, 0 on the boundary, -1 outside
 *   point_in_convex_polygon(poly, p): convex polygon in counter-clockwise order with n >= 3, O(log n), same return values
 *       collinear vertices allowed, as long as not all vertices are collinear
 *   simple_polygon(points): an order of the indices forming a simple polygon, O(n log n)
 *       points must be distinct and not all collinear
 *
 * Floating point part:
 *   rotate(center, p, degrees): rotates p counter-clockwise around center
 *   clockwise_angle(a, b): clockwise angle in degrees in [0, 360) from vector a to vector b
 *   great_circle_distance(lat1, lon1, lat2, lon2, radius): haversine formula, angles in degrees
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;

    Point operator-(const Point& p) const{
        return {x - p.x, y - p.y};
    }

    bool operator==(const Point& p) const{
        return x == p.x && y == p.y;
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

struct PointF{
    double x, y;
};

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
    Point a{0, 0}, b{4, 0}, c{4, 4}, d{0, 4};
    assert(cross(a, b, c) == 16 && cross(a, c, b) == -16 && cross(a, b, Point{8, 0}) == 0);
    assert(on_segment({2, 0}, a, b) && on_segment(a, a, b) && !on_segment({5, 0}, a, b) && !on_segment({2, 1}, a, b));
    assert(on_segment({3, 3}, {3, 3}, {3, 3}) && !on_segment({3, 4}, {3, 3}, {3, 3}));

    assert(segments_intersect(a, c, b, d) && segments_intersect(a, b, b, c) && !segments_intersect(a, b, c, d));
    assert(segments_intersect(a, b, {2, 0}, {6, 0}) && !segments_intersect(a, b, {5, 0}, {6, 0}));

    assert(abs(dist_point_segment({2, 3}, a, b) - 3.0) < 1e-9);
    assert(abs(dist_point_segment({7, 4}, a, b) - 5.0) < 1e-9);
    assert(abs(dist_point_segment({1, 1}, a, a) - sqrt(2.0)) < 1e-9);

    vector<Point> square = {a, b, c, d};
    assert(area2(square) == 32 && area2({d, c, b, a}) == -32);
    assert(point_in_convex_polygon(square, {2, 2}) == 1 && point_in_polygon(square, {2, 2}) == 1);
    assert(point_in_convex_polygon(square, {4, 2}) == 0 && point_in_polygon(square, {4, 2}) == 0);
    assert(point_in_convex_polygon(square, {0, 0}) == 0 && point_in_polygon(square, {0, 0}) == 0);
    assert(point_in_convex_polygon(square, {5, 2}) == -1 && point_in_polygon(square, {5, 2}) == -1);
    assert(point_in_convex_polygon(square, {0, 6}) == -1 && point_in_polygon(square, {0, 6}) == -1);

    vector<Point> notch = {{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}};
    assert(point_in_polygon(notch, {2, 3}) == -1 && point_in_polygon(notch, {2, 0}) == 0 && point_in_polygon(notch, {1, 1}) == 1);

    const long long E9 = 1000000000;
    vector<Point> huge = {{-E9, -E9}, {E9, -E9}, {E9, E9}, {-E9, E9}};
    assert(area2(huge) == 8 * E9 * E9);
    assert(point_in_convex_polygon(huge, {E9, E9}) == 0 && point_in_convex_polygon(huge, {0, 0}) == 1);

    vector<Point> pts = {{0, 0}, {2, 0}, {1, 1}, {0, 2}, {2, 2}};
    vector<int> order = simple_polygon(pts);
    assert((int)order.size() == 5 && order[0] == 0);

    PointF r = rotate({1, 1}, {2, 1}, 90);
    assert(abs(r.x - 1) < 1e-9 && abs(r.y - 2) < 1e-9);
    assert(abs(clockwise_angle({1, 0}, {0, -1}) - 90) < 1e-9 && abs(clockwise_angle({1, 0}, {0, 1}) - 270) < 1e-9);

    long double pi = acosl(-1.0L);
    assert(fabsl(great_circle_distance(0, 0, 0, 90, 1) - pi / 2) < 1e-12);
    assert(fabsl(great_circle_distance(0, 0, 0, 180, 1) - pi) < 1e-9);
    assert(fabsl(great_circle_distance(10, 20, 10, 20, 6371) - 0) < 1e-12);

    return 0;
}
