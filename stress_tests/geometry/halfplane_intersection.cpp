#include "../common.h"

#define main library_main
#include "../../code_library/geometry/halfplane_intersection.cpp"
#undef main

long double area2(const vector<Point>& poly){
    long double res = 0;
    for (int i = 0, n = poly.size(); i < n; i++) res += cross(poly[i], poly[(i + 1) % n]);
    return res;
}

/// Sutherland-Hodgman: the box square clipped by one half-plane at a time, O(n) per half-plane
vector<Point> brute_clip(const vector<Halfplane>& h, long double box){
    vector<Point> poly = {Point(-box, -box), Point(box, -box), Point(box, box), Point(-box, box)};

    for (auto& hp : h){
        vector<Point> next;
        for (int i = 0, n = poly.size(); i < n; i++){
            Point a = poly[i], b = poly[(i + 1) % n];
            long double sa = cross(hp.pq, a - hp.p), sb = cross(hp.pq, b - hp.p);
            if (sa >= 0) next.push_back(a);
            if ((sa > 0 && sb < 0) || (sa < 0 && sb > 0)) next.push_back(a + (b - a) * (sa / (sa - sb)));
        }
        poly = next;
    }
    return poly;
}

long double signed_dist(const Point& a, const Point& b, const Point& p){
    return cross(b - a, p - a) / hypotl(b.x - a.x, b.y - a.y);
}

/// Both polygons contain each other's vertices within tol, so they are the same region; the fast one must also be clean
void check(const vector<Halfplane>& h, long double box){
    vector<Point> fast = halfplane_intersection(h, box), slow = brute_clip(h, box);
    long double tol = 1e-6L, area_slow = fabsl(area2(slow)) / 2, perimeter = tol;
    for (int i = 0, n = slow.size(); i < n; i++) perimeter += hypotl(slow[i].x - slow[(i + 1) % n].x, slow[i].y - slow[(i + 1) % n].y);
    int k = fast.size();

    /// A zero-area region has rounding noise of about tol times its perimeter, plus an absolute floor for a single point inside a huge box
    if (k == 0){
        assert(area_slow <= tol * perimeter + 1e-9L);
        return;
    }

    assert(k >= 3 && area_slow > tol * perimeter);
    for (int i = 0; i < k; i++){
        const Point &a = fast[i], &b = fast[(i + 1) % k], &c = fast[(i + 2) % k];
        assert(hypotl(b.x - a.x, b.y - a.y) > EPS && signed_dist(a, b, c) > 0);
        assert(fabsl(a.x) <= box + tol && fabsl(a.y) <= box + tol);
        for (auto& hp : h) assert(cross(hp.pq, a - hp.p) / hp.len >= -tol);
        for (auto& p : slow) assert(signed_dist(a, b, p) >= -tol);
    }
    assert(fabsl(fabsl(area2(fast)) / 2 - area_slow) <= tol * perimeter);
}

Point random_point(int range){
    return Point(stress::rand_int(-range, range), stress::rand_int(-range, range));
}

/// A random line, oriented to keep the anchor when anchor is given, so most intersections stay non-empty
Halfplane random_halfplane(int range, const Point* anchor){
    Point a = random_point(range), b = random_point(range);
    while (a.x == b.x && a.y == b.y) b = random_point(range);
    if (anchor && cross(b - a, *anchor - a) < 0) swap(a, b);
    return Halfplane(a, b);
}

int main(){
    vector<Halfplane> family;
    int dirs[8][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
    for (int x = -1; x <= 1; x++) for (int y = -1; y <= 1; y++) for (auto& d : dirs) family.emplace_back(Point(x, y), Point(x + d[0], y + d[1]));

    for (auto& a : family) for (auto& b : family){
        check({a, b}, 3);
        check({a, b, family[stress::rand_int(0, family.size() - 1)]}, 2);
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, it % 10 ? 8 : 60), range = it % 3 ? 5 : vector<int>{1000, 10000}[stress::rand_int(0, 1)];
        long double box = vector<long double>{3, 50, 1e4, 1e8, 1e9}[stress::rand_int(0, 4)];
        Point anchor = random_point(range);
        bool feasible = stress::rand_int(0, 1);

        vector<Halfplane> h;
        for (int i = 0; i < n; i++){
            if (i && stress::rand_int(0, 5) == 0) h.push_back(h[stress::rand_int(0, i - 1)]);
            else h.push_back(random_halfplane(range, feasible ? &anchor : nullptr));
        }
        check(h, box);
    }

    /// Lines through one far point: either an unbounded cone cut by a huge box, or only that point, which must come out empty
    for (long long it = 0; it < stress::scaled(5000); it++){
        Point center = random_point(10000);
        vector<Halfplane> h;
        for (int i = 0, n = stress::rand_int(1, 6); i < n; i++){
            Point d = random_point(5);
            while (d.x == 0 && d.y == 0) d = random_point(5);
            h.emplace_back(center, center + d);
        }
        check(h, it % 2 ? 1e9 : 1e8);
    }

    for (long long it = 0; it < stress::scaled(30); it++){
        int n = stress::rand_int(500, 2000);
        Point anchor = random_point(100);
        vector<Halfplane> h;
        for (int i = 0; i < n; i++) h.push_back(random_halfplane(10000, &anchor));
        check(h, it % 2 ? 1e9 : 1e4);
    }

    /// 2e5 half-planes all containing [-1, 1]^2, plus its four sides: the answer is that square
    vector<Halfplane> big = {Halfplane(Point(-1, -1), Point(1, -1)), Halfplane(Point(1, -1), Point(1, 1)), Halfplane(Point(1, 1), Point(-1, 1)), Halfplane(Point(-1, 1), Point(-1, -1))};
    while (big.size() < 200000){
        Halfplane hp = random_halfplane(10000, nullptr);
        bool covers = true;
        for (int sx = -1; sx <= 1; sx += 2) for (int sy = -1; sy <= 1; sy += 2) covers &= cross(hp.pq, Point(sx, sy) - hp.p) >= 0;
        if (covers) big.push_back(hp);
    }
    shuffle(big.begin(), big.end(), stress::rng());
    vector<Point> square = halfplane_intersection(big, 1e9);
    assert(square.size() == 4 && fabsl(area2(square) - 8) < 1e-6);

    return 0;
}
