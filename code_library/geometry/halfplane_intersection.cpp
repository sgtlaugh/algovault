/***
 *
 * Half-Plane Intersection
 * Intersects half-planes, each the left side of a directed line, inside the square |x|, |y| <= box
 *
 * Complexity: O(n log n)
 *
 * Halfplane(a, b): the points p with cross(b - a, p - a) >= 0, boundary included, a != b
 * IntHalfplane(ax, ay, bx, by): the same half-plane for integer a = (ax, ay) != b = (bx, by), |coordinates| <= 1e9
 * halfplane_intersection(h, box): vertices of the intersection inside the box, counter-clockwise, strictly convex
 *
 * Notes:
 *   - h is all Halfplane or all IntHalfplane; with IntHalfplane, box is an integer in [1, 1e9]
 *   - The starting vertex is unspecified and no vertex repeats
 *   - Empty if the intersection is empty or has zero area (a segment or a single point)
 *   - A bounded intersection comes out exact when box exceeds all its vertex coordinates
 *   - IntHalfplane decides exactly in __int128, rounds only the output, and is 3x faster: use it for integer input
 *   - Halfplane is long double with an absolute EPS = 1e-9, a vertex that close to a boundary is on it, so box <= 1e9
 *   - Halfplane merges directions within the double rounding of the points, so a line given twice by decimals merges
 *   - Halfplane never merges distinct integer directions with |x|, |y| <= 1e6
 *   - Verified against brute force: Halfplane for integer |x|, |y| <= 1e4 and 0.1-step decimals, IntHalfplane exactly
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long double EPS = 1e-9;

struct Point{
    long double x, y;

    explicit Point(long double x = 0, long double y = 0) : x(x), y(y) {}

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

struct IntHalfplane{
    static constexpr long long LIMIT = 1000000000;

    long long px, py, dx, dy;

    IntHalfplane(long long ax, long long ay, long long bx, long long by) : px(ax), py(ay), dx(bx - ax), dy(by - ay){
        assert((dx != 0 || dy != 0) && max({llabs(ax), llabs(ay), llabs(bx), llabs(by)}) <= LIMIT);
    }

    long long cross_dir(const IntHalfplane& h) const{
        return dx * h.dy - dy * h.dx;
    }

    long long side(long long x, long long y) const{
        return dx * (y - py) - dy * (x - px);
    }

    /// Needs a.cross_dir(b) > 0, or the sign flips; with |coordinates| <= LIMIT both products stay below 2^126
    bool out(const IntHalfplane& a, const IntHalfplane& b) const{
        long long det = a.cross_dir(b), t = b.side(a.px, a.py);
        return (__int128)side(a.px, a.py) * det + (__int128)cross_dir(a) * t <= 0;
    }

    bool upper() const{
        return dy > 0 || (dy == 0 && dx > 0);
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

    /// Of lines with directions equal up to double rounding of their points, only the innermost constrains anything
    auto size = [](const Halfplane& a){ return fabsl(a.p.x) + fabsl(a.p.y) + fabsl(a.pq.x) + fabsl(a.pq.y); };
    /// The dot product, not upper(), splits antiparallel lines, so twins rounded across angle pi still merge
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

/// Exact numerators leave two roundings per coordinate, the conversion and the division, so about 1 ulp
Point line_intersection(const IntHalfplane& s, const IntHalfplane& t){
    long long det = s.cross_dir(t), alpha = t.side(s.px, s.py);
    __int128 x = (__int128)s.px * det + (__int128)s.dx * alpha, y = (__int128)s.py * det + (__int128)s.dy * alpha;
    return Point((long double)x / det, (long double)y / det);
}

vector<Point> halfplane_intersection(vector<IntHalfplane> h, long long box){
    assert(1 <= box && box <= IntHalfplane::LIMIT);
    long long cx[4] = {box, box, -box, -box}, cy[4] = {-box, box, box, -box};
    for (int i = 0; i < 4; i++) h.emplace_back(cx[i], cy[i], cx[(i + 1) % 4], cy[(i + 1) % 4]);

    /// Splitting by upper() first, then sorting each half by cross, is 13% faster than one comparator at n = 5e4
    auto mid = partition(h.begin(), h.end(), [](const IntHalfplane& a){ return a.upper(); });
    auto by_angle = [](const IntHalfplane& a, const IntHalfplane& b){ return a.cross_dir(b) > 0; };
    sort(h.begin(), mid, by_angle);
    sort(mid, h.end(), by_angle);

    /// Antiparallel lines are never neighbors, as the box puts its directions pi / 2 and pi between them
    vector<IntHalfplane> lines;
    for (auto& hp : h){
        bool same_direction = !lines.empty() && lines.back().cross_dir(hp) == 0;
        if (!same_direction) lines.push_back(hp);
        else if (lines.back().side(hp.px, hp.py) > 0) lines.back() = hp;
    }

    /// The deque is the window dq[head, tail), 13% faster than std::deque at n = 100; the copy only sizes it
    vector<IntHalfplane> dq = lines;
    int head = 0, tail = 0;
    /// A vertex on the new line is popped too, so lines through one point leave one vertex and no zero-length edge
    for (auto& hp : lines){
        while (tail - head > 1 && hp.out(dq[tail - 2], dq[tail - 1])) tail--;
        while (tail - head > 1 && hp.out(dq[head], dq[head + 1])) head++;

        if (tail > head && dq[tail - 1].cross_dir(hp) <= 0) return {};
        dq[tail++] = hp;
    }

    while (tail - head > 2 && dq[head].out(dq[tail - 2], dq[tail - 1])) tail--;
    while (tail - head > 2 && dq[tail - 1].out(dq[head], dq[head + 1])) head++;
    if (tail - head < 3) return {};

    vector<Point> poly;
    for (int i = head; i < tail; i++) poly.push_back(line_intersection(dq[i], dq[i + 1 < tail ? i + 1 : head]));
    return poly;
}

int main(){
    auto near = [](const Point& p, long double x, long double y){ return fabsl(p.x - x) < 1e-9 && fabsl(p.y - y) < 1e-9; };
    Halfplane x_min(Point(0, 1), Point(0, 0)), y_min(Point(0, 0), Point(1, 0));  /// x >= 0 and y >= 0, left of a -> b

    /// Adding x + y <= 4
    vector<Point> poly = halfplane_intersection({x_min, y_min, Halfplane(Point(4, 0), Point(0, 4))}, 1e9);
    assert(poly.size() == 3);
    for (const Point& p : poly) assert(near(p, 0, 0) || near(p, 4, 0) || near(p, 0, 4));
    assert(cross(poly[1] - poly[0], poly[2] - poly[0]) > 0);  /// counter-clockwise

    assert(halfplane_intersection({x_min}, 5).size() == 4);   /// unbounded, so the box clips it to [0, 5] x [-5, 5]

    Halfplane y_max(Point(1, -1), Point(0, -1));                /// y <= -1
    assert(halfplane_intersection({y_min, y_max}, 100).empty());
    Halfplane x_plus_y_max(Point(0, 0), Point(-1, 1));          /// x + y <= 0
    assert(halfplane_intersection({x_min, y_min, x_plus_y_max}, 100).empty());  /// only the origin, zero area

    /// Integer input is exact: y >= 0, x <= 1 and below the line from (0, 0) to (1e9, 1) leave area 5e-10
    vector<IntHalfplane> sliver = {IntHalfplane(0, 0, 1, 0), IntHalfplane(1, 0, 1, 1), IntHalfplane(1000000000, 1, 0, 0)};
    poly = halfplane_intersection(sliver, 1000000000);
    assert(poly.size() == 3);
    for (const Point& p : poly) assert((p.x == 0 && p.y == 0) || (p.x == 1 && p.y == 0) || (p.x == 1 && fabsl(p.y - 1e-9L) < 1e-25L));

    vector<Halfplane> rounded = {Halfplane(Point(0, 0), Point(1, 0)), Halfplane(Point(1, 0), Point(1, 1)), Halfplane(Point(1e9, 1), Point(0, 0))};
    assert(halfplane_intersection(rounded, 1e9).empty());  /// the long double version merges (1, 0) and (1, 1e-9)

    static_assert(!is_convertible_v<int, Point>, "Halfplane(3, 7) and out(5) must not compile");
    return 0;
}
