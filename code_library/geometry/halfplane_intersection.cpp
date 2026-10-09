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
 *     Directions are sorted and compared for parallelism exactly, so integer inputs never mix up parallel lines
 *     Verified against brute force clipping for integer coordinates |x|, |y| <= 1e4 with box up to 1e9
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
    vector<Halfplane> lines;
    for (auto& hp : h){
        if (!lines.empty() && lines.back().upper() == hp.upper() && cross(lines.back().pq, hp.pq) == 0){
            if (cross(lines.back().pq, hp.p - lines.back().p) > 0) lines.back() = hp;
            continue;
        }
        lines.push_back(hp);
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
    auto matches = [](const vector<Point>& got, const vector<Point>& expected){
        int n = expected.size();
        if ((int)got.size() != n) return false;

        for (int shift = 0; shift < n; shift++){
            bool ok = true;
            for (int i = 0; i < n && ok; i++){
                const Point &a = got[(i + shift) % n], &b = expected[i];
                ok = fabsl(a.x - b.x) < 1e-6 && fabsl(a.y - b.y) < 1e-6;
            }
            if (ok) return true;
        }
        return false;
    };

    Halfplane x_min(Point(0, 1), Point(0, 0)), y_min(Point(0, 0), Point(1, 0));

    assert(matches(halfplane_intersection({x_min, y_min, Halfplane(Point(4, 0), Point(0, 4))}, 1e9), {Point(0, 0), Point(4, 0), Point(0, 4)}));
    assert(matches(halfplane_intersection({}, 10), {Point(-10, -10), Point(10, -10), Point(10, 10), Point(-10, 10)}));
    assert(matches(halfplane_intersection({x_min}, 5), {Point(0, -5), Point(5, -5), Point(5, 5), Point(0, 5)}));

    assert(matches(halfplane_intersection({y_min, Halfplane(Point(0, -1), Point(1, -1))}, 100), {Point(-100, 0), Point(100, 0), Point(100, 100), Point(-100, 100)}));
    assert(halfplane_intersection({Halfplane(Point(0, 1), Point(1, 1)), Halfplane(Point(1, 0), Point(0, 0))}, 100).empty());
    assert(halfplane_intersection({y_min, Halfplane(Point(1, 0), Point(0, 0))}, 100).empty());
    assert(halfplane_intersection({x_min, y_min, Halfplane(Point(0, 0), Point(-1, 1))}, 100).empty());

    vector<Halfplane> square = {x_min, y_min, Halfplane(Point(2, 0), Point(2, 1)), Halfplane(Point(1, 2), Point(0, 2))};
    vector<Halfplane> crowded = square;
    crowded.insert(crowded.end(), square.begin(), square.end());
    crowded.push_back(Halfplane(Point(3, 0), Point(3, 1)));
    crowded.push_back(Halfplane(Point(0, 0), Point(1, -1)));
    crowded.push_back(Halfplane(Point(2, 2), Point(1, 3)));
    assert(matches(halfplane_intersection(crowded, 1e9), {Point(0, 0), Point(2, 0), Point(2, 2), Point(0, 2)}));

    vector<Halfplane> cut = {Halfplane(Point(2, 0), Point(3, 1)), Halfplane(Point(0, 4), Point(0, 0)), Halfplane(Point(6, 0), Point(0, 6)),
                             Halfplane(Point(4, 4), Point(0, 4)), Halfplane(Point(0, 0), Point(4, 0)), Halfplane(Point(4, 0), Point(4, 4))};
    assert(matches(halfplane_intersection(cut, 1e9), {Point(0, 0), Point(2, 0), Point(4, 2), Point(2, 4), Point(0, 4)}));

    return 0;
}
