/***
 *
 * Half-plane Intersection
 * Intersects half-planes, each the left side of a directed line, inside the square |x|, |y| <= box
 *
 * Complexity: O(n log n)
 *
 * Halfplane(a, b) keeps the points p with cross(b - a, p - a) >= 0, boundary included, a != b
 * halfplane_intersection(h, box) returns the vertices of the intersection clipped to the box, counter-clockwise,
 *     starting vertex unspecified, strictly convex and without repeated vertices
 *     Empty if the intersection is empty or has zero area (a segment or a single point)
 *     A bounded intersection comes out exact when box exceeds all its vertex coordinates
 *
 * Precision: long double, EPS = 1e-9 is a distance tolerance, a vertex within EPS of a boundary counts as on it
 *     Directions are sorted exactly, and two count as one when they differ by no more than the double rounding of the points,
 *     so a line given twice through different decimal points is merged, while distinct integer directions with |x|, |y| <= 1e6 never are
 *     Verified against brute force clipping for integer coordinates |x|, |y| <= 1e4 with box up to 1e9, and for 0.1-step decimals
 *     EPS is absolute, about 9 long double ulps at box 1e9, so a box much beyond 1e9 is not supported
 *
 * Example, x >= 0, y >= 0, x + y <= 4:
 *     vector<Halfplane> h = {Halfplane(Point(0, 0), Point(1, 0)), Halfplane(Point(0, 1), Point(0, 0)), Halfplane(Point(4, 0), Point(0, 4))};
 *     auto poly = halfplane_intersection(h, 1e9);  /// (0, 0), (4, 0), (0, 4)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long double EPS = 1e-9;

struct Point{
    long double x, y;

    Point(long double x = 0, long double y = 0) : x(x), y(y) {}

    Point operator+(const Point& p) const{
        return Point(x + p.x, y + p.y);
    }

    Point operator-(const Point& p) const{
        return Point(x - p.x, y - p.y);
    }

    Point operator*(long double k) const{
        return Point(x * k, y * k);
    }
};

long double cross(const Point& a, const Point& b){
    return a.x * b.y - a.y * b.x;
}

struct Halfplane{
    Point p, pq;
    long double len;

    Halfplane(const Point& a, const Point& b) : p(a), pq(b - a), len(hypotl(pq.x, pq.y)) {}

    bool out(const Point& r) const{
        return cross(pq, r - p) < -EPS * len;
    }

    /// Angle in [0, pi), tested on the exact direction so -0.0 and parallel lines stay consistent
    bool upper() const{
        return pq.y > 0 || (pq.y == 0 && pq.x > 0);
    }
};

Point line_intersection(const Halfplane& s, const Halfplane& t){
    long double alpha = cross(t.p - s.p, t.pq) / cross(s.pq, t.pq);
    return s.p + s.pq * alpha;
}

vector<Point> halfplane_intersection(vector<Halfplane> h, long double box){
    Point corners[4] = {Point(box, -box), Point(box, box), Point(-box, box), Point(-box, -box)};
    for (int i = 0; i < 4; i++) h.emplace_back(corners[i], corners[(i + 1) % 4]);

    sort(h.begin(), h.end(), [](const Halfplane& a, const Halfplane& b){
        if (a.upper() != b.upper()) return a.upper();
        return cross(a.pq, b.pq) > 0;
    });

    /// Of each group of same-direction lines only the innermost constrains anything
    /// One line given by two different decimal point pairs gets directions that differ by the double rounding of its points,
    /// a few DBL_EPSILON times the coordinate sizes, while distinct integer directions have cross >= 1, well above that up to 1e6
    /// The dot product, not upper(), separates antiparallel lines, so twins rounded to either side of angle pi still merge
    auto size = [](const Halfplane& a){ return fabsl(a.p.x) + fabsl(a.p.y) + fabsl(a.pq.x) + fabsl(a.pq.y); };
    auto same_direction = [&](const Halfplane& a, const Halfplane& b){
        return a.pq.x * b.pq.x + a.pq.y * b.pq.y > 0 && fabsl(cross(a.pq, b.pq)) <= 8 * DBL_EPSILON * (size(a) * b.len + size(b) * a.len);
    };
    auto merge = [](Halfplane& kept, const Halfplane& hp){
        if (cross(kept.pq, hp.p - kept.p) > 0) kept = hp;
    };
    vector<Halfplane> lines;
    for (auto& hp : h){
        if (!lines.empty() && same_direction(lines.back(), hp)) merge(lines.back(), hp);
        else lines.push_back(hp);
    }

    /// A direction rounded to just below angle 0 sorts last, next to its twin at the front
    if (same_direction(lines.back(), lines[0])){
        merge(lines[0], lines.back());
        lines.pop_back();
    }

    deque<Halfplane> dq;
    for (auto& hp : lines){
        while (dq.size() > 1 && hp.out(line_intersection(dq[dq.size() - 1], dq[dq.size() - 2]))) dq.pop_back();
        while (dq.size() > 1 && hp.out(line_intersection(dq[0], dq[1]))) dq.pop_front();

        /// Everything between two opposite directions was cut away, so nothing is left with positive area
        if (!dq.empty() && cross(dq.back().pq, hp.pq) <= 0) return {};
        dq.push_back(hp);
    }

    while (dq.size() > 2 && dq[0].out(line_intersection(dq[dq.size() - 1], dq[dq.size() - 2]))) dq.pop_back();
    while (dq.size() > 2 && dq.back().out(line_intersection(dq[0], dq[1]))) dq.pop_front();
    if (dq.size() < 3) return {};

    /// Lines through a common vertex leave zero-length edges, a degenerate intersection collapses to under 3 vertices
    auto same = [](const Point& a, const Point& b){ return fabsl(a.x - b.x) <= EPS && fabsl(a.y - b.y) <= EPS; };
    vector<Point> poly;
    for (int i = 0, k = dq.size(); i < k; i++){
        Point v = line_intersection(dq[i], dq[(i + 1) % k]);
        if (poly.empty() || !same(poly.back(), v)) poly.push_back(v);
    }
    while (poly.size() > 1 && same(poly.back(), poly[0])) poly.pop_back();

    if (poly.size() < 3) return {};
    return poly;
}

int main(){
    auto near = [](const Point& p, long double x, long double y){ return fabsl(p.x - x) < 1e-9 && fabsl(p.y - y) < 1e-9; };
    Halfplane x_min(Point(0, 1), Point(0, 0)), y_min(Point(0, 0), Point(1, 0));  /// x >= 0 and y >= 0, left of a -> b

    vector<Point> poly = halfplane_intersection({x_min, y_min, Halfplane(Point(4, 0), Point(0, 4))}, 1e9);  /// and x + y <= 4
    assert(poly.size() == 3);
    for (const Point& p : poly) assert(near(p, 0, 0) || near(p, 4, 0) || near(p, 0, 4));
    assert(cross(poly[1] - poly[0], poly[2] - poly[0]) > 0);  /// counter-clockwise

    assert(halfplane_intersection({x_min}, 5).size() == 4);   /// unbounded, so the box clips it to [0, 5] x [-5, 5]

    Halfplane y_max(Point(1, -1), Point(0, -1));                /// y <= -1
    assert(halfplane_intersection({y_min, y_max}, 100).empty());
    Halfplane x_plus_y_max(Point(0, 0), Point(-1, 1));          /// x + y <= 0
    assert(halfplane_intersection({x_min, y_min, x_plus_y_max}, 100).empty());  /// only the origin, zero area
    return 0;
}
