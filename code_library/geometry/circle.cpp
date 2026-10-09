/***
 *
 * Circle Geometry
 * Intersections, tangents, enclosing circles, covering and areas of circles in long double
 *
 * Complexity: O(1) for intersections, tangents, circumcircle and circle_circle_area
 *             minimum_enclosing_circle expected O(n), circle_polygon_area O(n)
 *             circle_union_area and max_circle_cover O(n^2 log n)
 *
 * Circle{center, r} with r >= 0, tangent and touching cases are detected with a relative tolerance EPS
 *   circle_line_intersection(c, a, b): where the line through a != b meets c, 0, 1 (tangent) or 2 points in the a -> b direction
 *       the tolerance scales with r, so r = 0 meets the line only when the center lies exactly on it
 *   circle_circle_intersection(a, b): 0, 1 (tangent) or 2 points, none for concentric circles (identical ones included)
 *   tangents_from_point(c, p): tangency points of the lines through p touching c (r > 0), none if p is inside, p itself if on c
 *   common_tangents(a, b): every common tangent line as its tangency points {on a, on b}, r > 0 for both, outer ones first
 *       0 nested, 1 internally tangent, 2 overlapping, 3 externally tangent, 4 apart, none for identical circles
 *       when both points coincide the line goes through them perpendicular to the line of centers
 *   circumcircle(a, b, c): the circle through three non-collinear points
 *   minimum_enclosing_circle(points): smallest circle containing every point, points non-empty, shuffled Welzl
 *   circle_polygon_area(c, poly): area of c intersected with a simple polygon in either orientation
 *   circle_circle_area(a, b): area of the intersection of two circles
 *   circle_union_area(circles): area of the union, Green's theorem over the uncovered boundary arcs
 *   max_circle_cover(points, r): {most points one circle of radius r covers, a center achieving it}, boundary points count
 *
 * Verified with |coordinates|, radii <= 1e4
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

using Float = long double;

const Float EPS = 1e-12L;
const Float PI = acosl(-1.0L);

struct Point{
    Float x, y;

    Point operator+(const Point& p) const{
        return {x + p.x, y + p.y};
    }

    Point operator-(const Point& p) const{
        return {x - p.x, y - p.y};
    }

    Point operator*(Float k) const{
        return {x * k, y * k};
    }

    Point operator/(Float k) const{
        return {x / k, y / k};
    }

    Float dot(const Point& p) const{
        return x * p.x + y * p.y;
    }

    Float cross(const Point& p) const{
        return x * p.y - y * p.x;
    }

    Float norm2() const{
        return x * x + y * y;
    }

    Float norm() const{
        return sqrtl(norm2());
    }

    Point perp() const{
        return {-y, x};
    }
};

struct Circle{
    Point center;
    Float r;
};

vector<Point> circle_line_intersection(const Circle& c, const Point& a, const Point& b){
    Point ab = b - a;
    Float len2 = ab.norm2(), s = ab.cross(c.center - a);
    Float h2 = c.r * c.r * len2 - s * s, tol = EPS * c.r * c.r * len2;
    if (h2 < -tol) return {};

    Point foot = a + ab * (ab.dot(c.center - a) / len2);
    if (h2 <= tol) return {foot};

    Point h = ab * (sqrtl(h2) / len2);
    return {foot - h, foot + h};
}

vector<Point> circle_circle_intersection(const Circle& a, const Circle& b){
    Point d = b.center - a.center;
    Float d2 = d.norm2(), sum = a.r + b.r, dif = a.r - b.r, tol = EPS * sum * sum;
    if (d2 == 0 || d2 > sum * sum + tol || d2 < dif * dif - tol) return {};

    Float p = (d2 + a.r * a.r - b.r * b.r) / (2 * d2);
    Point mid = a.center + d * p;
    if (fabsl(d2 - sum * sum) <= tol || fabsl(d2 - dif * dif) <= tol) return {mid};

    Point h = d.perp() * sqrtl(max(a.r * a.r - p * p * d2, (Float)0) / d2);
    return {mid + h, mid - h};
}

vector<Point> tangents_from_point(const Circle& c, const Point& p){
    Point d = p - c.center;
    Float d2 = d.norm2(), r2 = c.r * c.r;
    if (d2 < r2 * (1 - EPS)) return {};
    if (d2 <= r2 * (1 + EPS)) return {p};

    Point foot = c.center + d * (r2 / d2), h = d.perp() * (c.r * sqrtl(d2 - r2) / d2);
    return {foot + h, foot - h};
}

vector<pair<Point, Point>> common_tangents(const Circle& a, const Circle& b){
    vector<pair<Point, Point>> res;
    Point d = b.center - a.center;
    Float d2 = d.norm2();
    if (d2 == 0) return res;

    /// A negated radius puts the tangency point on the far side of b, turning outer tangents into inner ones
    for (Float rb : {b.r, -b.r}){
        Float dr = a.r - rb, h2 = d2 - dr * dr;
        if (h2 < -EPS * d2) continue;

        Float h = h2 <= EPS * d2 ? 0 : sqrtl(h2);
        for (int sign : {1, -1}){
            Point v = (d * dr + d.perp() * (h * sign)) / d2;
            res.push_back({a.center + v * a.r, b.center + v * rb});
            if (h == 0) break;
        }
    }

    return res;
}

Circle circumcircle(const Point& a, const Point& b, const Point& c){
    Point u = b - a, v = c - a;
    Point center = a + (v * u.norm2() - u * v.norm2()).perp() / (v.cross(u) * 2);
    return {center, max({(center - a).norm(), (center - b).norm(), (center - c).norm()})};
}

Circle minimum_enclosing_circle(vector<Point> points){
    assert(!points.empty());
    shuffle(points.begin(), points.end(), mt19937(chrono::steady_clock::now().time_since_epoch().count()));

    int n = points.size();
    Circle c{points[0], 0};
    auto outside = [&](const Point& p){ return (p - c.center).norm() > c.r * (1 + EPS); };
    for (int i = 1; i < n; i++){
        if (!outside(points[i])) continue;
        c = {points[i], 0};
        for (int j = 0; j < i; j++){
            if (!outside(points[j])) continue;
            /// The rounded center can sit farther than |pi - pj| / 2 from pi or pj, and with a tiny radius that error beats the
            /// relative tolerance, so pi itself would test outside and reach circumcircle as a degenerate triple giving NaN
            Point mid = (points[i] + points[j]) / 2;
            c = {mid, max((mid - points[i]).norm(), (mid - points[j]).norm())};
            for (int k = 0; k < j; k++){
                if (outside(points[k])) c = circumcircle(points[i], points[j], points[k]);
            }
        }
    }

    return c;
}

Float circle_polygon_area(const Circle& c, const vector<Point>& poly){
    Float half_r2 = c.r * c.r / 2, res = 0;
    auto angle = [](const Point& p, const Point& q){ return atan2l(p.cross(q), p.dot(q)); };

    /// Signed area of the circle intersected with the triangle (center, p, q), center at the origin
    auto sector = [&](const Point& p, const Point& q) -> Float{
        Point d = q - p;
        Float len2 = d.norm2();
        if (len2 == 0) return 0;

        Float a = d.dot(p) / len2, b = (p.norm2() - c.r * c.r) / len2, det = a * a - b;
        if (det <= 0) return angle(p, q) * half_r2;

        Float s = max((Float)0, -a - sqrtl(det)), t = min((Float)1, -a + sqrtl(det));
        if (t < 0 || s >= 1) return angle(p, q) * half_r2;

        Point u = p + d * s, v = p + d * t;
        return angle(p, u) * half_r2 + u.cross(v) / 2 + angle(v, q) * half_r2;
    };

    for (int i = 0, n = poly.size(); i < n; i++) res += sector(poly[i] - c.center, poly[(i + 1) % n] - c.center);
    return fabsl(res);
}

Float circle_circle_area(const Circle& a, const Circle& b){
    Float d = (a.center - b.center).norm(), r = min(a.r, b.r);
    if (d >= a.r + b.r) return 0;
    if (d <= fabsl(a.r - b.r)) return PI * r * r;

    Float alpha = acosl(clamp((d * d + a.r * a.r - b.r * b.r) / (2 * d * a.r), (Float)-1, (Float)1));
    Float beta = acosl(clamp((d * d + b.r * b.r - a.r * a.r) / (2 * d * b.r), (Float)-1, (Float)1));
    return a.r * a.r * (alpha - sinl(2 * alpha) / 2) + b.r * b.r * (beta - sinl(2 * beta) / 2);
}

Float circle_union_area(const vector<Circle>& circles){
    int n = circles.size();
    Float res = 0;
    for (int i = 0; i < n; i++){
        const Circle& c = circles[i];
        bool hidden = c.r <= 0;
        vector<pair<Float, Float>> covered;
        for (int j = 0; j < n && !hidden; j++){
            const Circle& o = circles[j];
            Point v = o.center - c.center;
            Float d = v.norm();

            /// Of identical circles only the first one counts, otherwise each would hide the other
            if (j != i && d <= o.r - c.r && (o.r > c.r || j < i)) hidden = true;
            if (j == i || d >= c.r + o.r || d <= fabsl(c.r - o.r)) continue;

            Float mid = atan2l(v.y, v.x), half = acosl(clamp((c.r * c.r + d * d - o.r * o.r) / (2 * c.r * d), (Float)-1, (Float)1));
            Float lo = mid - half, hi = mid + half;
            if (lo < 0) lo += 2 * PI, hi += 2 * PI;
            if (hi <= 2 * PI) covered.push_back({lo, hi});
            else covered.push_back({lo, 2 * PI}), covered.push_back({0, hi - 2 * PI});
        }
        if (hidden) continue;

        /// Green's theorem: the arc from angle s to t adds half the integral of x dy - y dx along it
        auto arc = [&](Float s, Float t){
            return (c.r * c.r * (t - s) + c.center.x * c.r * (sinl(t) - sinl(s)) - c.center.y * c.r * (cosl(t) - cosl(s))) / 2;
        };

        sort(covered.begin(), covered.end());
        Float from = 0;
        for (auto [lo, hi] : covered){
            if (lo > from) res += arc(from, lo);
            from = max(from, hi);
        }
        if (from < 2 * PI) res += arc(from, 2 * PI);
    }

    return res;
}

pair<int, Point> max_circle_cover(const vector<Point>& points, Float r){
    int n = points.size(), best = 0;
    Point center{0, 0};
    for (int i = 0; i < n; i++){
        int count = 0;
        vector<pair<Float, int>> events;
        for (int j = 0; j < n; j++){
            Point v = points[j] - points[i];
            Float d2 = v.norm2();
            if (d2 == 0){
                count++;
                continue;
            }
            if (d2 > 4 * r * r * (1 + EPS)) continue;

            /// Centers on the circle of radius r around points[i] that also cover points[j], widened by EPS for boundary ties
            Float mid = atan2l(v.y, v.x), half = atan2l(sqrtl(max((Float)0, 4 * r * r - d2)), sqrtl(d2));
            Float lo = mid - half - EPS, hi = mid + half + EPS;
            if (lo < 0) lo += 2 * PI, hi += 2 * PI;
            events.push_back({lo, -1});
            events.push_back({min(hi, 2 * PI), 1});
            if (hi > 2 * PI) events.push_back({0, -1}), events.push_back({hi - 2 * PI, 1});
        }

        /// -1 opens an interval, so at equal angles it sorts before the closing +1
        sort(events.begin(), events.end());
        Float at = 0;
        int most = count;
        for (auto [angle, type] : events){
            count -= type;
            if (count > most) most = count, at = angle;
        }
        if (most > best) best = most, center = points[i] + Point{cosl(at), sinl(at)} * r;
    }

    return {best, center};
}

int main(){
    auto near = [](const Point& p, Float x, Float y){ return fabsl(p.x - x) < 1e-9 && fabsl(p.y - y) < 1e-9; };

    vector<Point> pts = circle_line_intersection({{0, 0}, 5}, {-10, 3}, {10, 3});
    assert(near(pts[0], -4, 3) && near(pts[1], 4, 3));         /// ordered from a towards b
    assert(circle_line_intersection({{0, 0}, 5}, {10, 5}, {-10, 5}).size() == 1);  /// tangent at (0, 5)

    pts = circle_circle_intersection({{0, 0}, 5}, {{8, 0}, 5});
    assert(near(pts[0], 4, 3) && near(pts[1], 4, -3));

    pts = tangents_from_point({{0, 0}, 3}, {5, 0});
    assert(near(pts[0], 1.8, 2.4) && near(pts[1], 1.8, -2.4));
    assert(tangents_from_point({{0, 0}, 3}, {1, 1}).empty());  /// p is inside

    auto tangents = common_tangents({{0, 0}, 1}, {{4, 0}, 1});
    assert(tangents.size() == 4);                              /// apart: 2 outer and 2 inner
    assert(near(tangents[0].first, 0, 1) && near(tangents[0].second, 4, 1));  /// outer ones first

    Circle c = circumcircle({0, 0}, {4, 0}, {0, 3});
    assert(near(c.center, 2, 1.5) && fabsl(c.r - 2.5) < 1e-9);  /// the hypotenuse is a diameter

    c = minimum_enclosing_circle({{0, 0}, {2, 0}, {0, 2}, {2, 2}, {1, 1}});
    assert(near(c.center, 1, 1) && fabsl(c.r - sqrtl(2)) < 1e-9);

    vector<Point> square = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
    assert(fabsl(circle_polygon_area({{0, 0}, 1}, square) - PI) < 1e-9);  /// the circle fits inside

    Float lens = 2 * PI / 3 - sqrtl(3) / 2;                    /// unit circles one radius apart
    assert(fabsl(circle_circle_area({{0, 0}, 1}, {{1, 0}, 1}) - lens) < 1e-9);
    assert(fabsl(circle_union_area({{{0, 0}, 1}, {{1, 0}, 1}}) - (2 * PI - lens)) < 1e-9);

    auto [covered, at] = max_circle_cover({{0, 0}, {2, 0}, {4, 0}, {10, 10}}, 2);
    assert(covered == 3);
    assert(near(at, 2, 0));
    return 0;
}
