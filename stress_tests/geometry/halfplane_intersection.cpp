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

    /// Zero-area rounding noise is about tol times the perimeter, plus a floor for a single point inside a huge box
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

__int128 gcd128(__int128 a, __int128 b){
    while (b){
        __int128 r = a % b;
        a = b;
        b = r;
    }

    return a < 0 ? -a : a;
}

struct Corner{
    __int128 x, y, d;
    vector<int> tight;
};

/// Exact rationals: a corner is a line crossing inside every half-plane where two non-parallel lines are tight
vector<Corner> brute_corners(const vector<IntHalfplane>& h){
    vector<Corner> res;
    for (int i = 0, n = h.size(); i < n; i++) for (int j = i + 1; j < n; j++){
        const IntHalfplane &a = h[i], &b = h[j];
        __int128 d = (__int128)a.dx * b.dy - (__int128)a.dy * b.dx;
        if (d == 0) continue;

        __int128 t = (__int128)(b.px - a.px) * b.dy - (__int128)(b.py - a.py) * b.dx;
        __int128 x = a.px * d + a.dx * t, y = a.py * d + a.dy * t;
        if (d < 0){
            x = -x, y = -y, d = -d;
        }
        __int128 g = gcd128(gcd128(x, y), d);
        Corner c = {x / g, y / g, d / g, {}};

        bool inside = true;
        for (int k = 0; k < n && inside; k++){
            __int128 s = h[k].dx * (c.y - h[k].py * c.d) - h[k].dy * (c.x - h[k].px * c.d);
            inside = s >= 0;
            if (s == 0) c.tight.push_back(k);
        }
        bool turns = false;
        for (int k : c.tight) turns |= h[k].cross_dir(h[c.tight[0]]) != 0;
        bool seen = false;
        for (auto& o : res) seen |= o.x == c.x && o.y == c.y && o.d == c.d;
        if (inside && turns && !seen) res.push_back(c);
    }

    return res;
}

/// Same corner set exactly as the brute force, each output vertex rounded from its corner, in boundary order
void check_exact(const vector<IntHalfplane>& h, long long box){
    vector<IntHalfplane> all = h;
    long long cx[4] = {box, box, -box, -box}, cy[4] = {-box, box, box, -box};
    for (int i = 0; i < 4; i++) all.emplace_back(cx[i], cy[i], cx[(i + 1) % 4], cy[(i + 1) % 4]);

    vector<Point> fast = halfplane_intersection(h, box);
    vector<Corner> corners = brute_corners(all);
    int k = fast.size();

    if (corners.size() < 3){
        assert(k == 0);
        return;
    }

    assert(k == (int)corners.size());
    vector<int> match(k);
    vector<bool> used(k);
    for (int i = 0; i < k; i++){
        long double best = INFINITY;
        for (int j = 0; j < k; j++){
            long double dist = hypotl(fast[i].x - (long double)corners[j].x / corners[j].d, fast[i].y - (long double)corners[j].y / corners[j].d);
            if (dist < best){
                best = dist;
                match[i] = j;
            }
        }
        assert(best <= 1e-15L * (box + 1) && !used[match[i]]);
        used[match[i]] = true;
    }

    /// Consecutive vertices share a tight line, so the output walks the boundary instead of jumping across it
    vector<int> edge(k, -1);
    for (int i = 0; i < k; i++){
        const vector<int> &a = corners[match[i]].tight, &b = corners[match[(i + 1) % k]].tight;
        for (int x : a) if (find(b.begin(), b.end(), x) != b.end()) edge[i] = x;
        assert(edge[i] != -1);
    }

    /// Exact counter-clockwise test: each edge turns left from the one before, where a clockwise walk turns right
    for (int i = 0; i < k; i++) assert(all[edge[(i + k - 1) % k]].cross_dir(all[edge[i]]) > 0);
}

IntHalfplane random_int_halfplane(long long range, const pair<long long, long long>* anchor){
    long long ax = stress::rand_int(-range, range), ay = stress::rand_int(-range, range), bx = ax, by = ay;
    while (bx == ax && by == ay){
        bx = stress::rand_int(-range, range);
        by = stress::rand_int(-range, range);
    }
    IntHalfplane hp(ax, ay, bx, by);
    if (anchor && hp.side(anchor->first, anchor->second) < 0) hp = IntHalfplane(bx, by, ax, ay);
    return hp;
}

void stress_exact(){
    vector<IntHalfplane> family;
    int dirs[8][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
    for (int x = -1; x <= 1; x++) for (int y = -1; y <= 1; y++) for (auto& d : dirs) family.emplace_back(x, y, x + d[0], y + d[1]);

    for (auto& a : family) for (auto& b : family){
        check_exact({a, b}, 1);
        check_exact({a, b, family[stress::rand_int(0, family.size() - 1)]}, 2);
    }

    /// Tiny ranges force parallel twins, antiparallel strips, lines through one point and zero-area answers
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, it % 10 ? 8 : 40);
        long long range = vector<long long>{2, 5, 1000, 1000000000}[it % 4];
        long long box = it % 4 == 3 ? IntHalfplane::LIMIT - stress::rand_int(0, 1) : vector<long long>{1, 3, 50, 1000}[stress::rand_int(0, 3)];
        pair<long long, long long> anchor = {stress::rand_int(-range, range), stress::rand_int(-range, range)};
        bool feasible = stress::rand_int(0, 1);

        vector<IntHalfplane> h;
        for (int i = 0; i < n; i++){
            if (i && stress::rand_int(0, 5) == 0) h.push_back(h[stress::rand_int(0, i - 1)]);
            else h.push_back(random_int_halfplane(range, feasible ? &anchor : nullptr));
        }
        check_exact(h, box);
    }

    /// Lines through one point at the coordinate limit: a cone, or only that point, which has zero area
    for (long long it = 0; it < stress::scaled(5000); it++){
        long long lim = IntHalfplane::LIMIT, x = stress::rand_int(-lim, lim), y = stress::rand_int(-lim, lim);
        if (it % 2) x = stress::rand_int(0, 1) ? lim : -lim;
        vector<IntHalfplane> h;
        for (int i = 0, n = stress::rand_int(1, 6); i < n; i++){
            long long dx = 0, dy = 0;
            while ((dx == 0 && dy == 0) || llabs(x + dx) > lim || llabs(y + dy) > lim){
                dx = stress::rand_int(-3, 3);
                dy = stress::rand_int(-3, 3);
            }
            h.emplace_back(x, y, x + dx, y + dy);
        }
        check_exact(h, lim);
    }

    /// Triangles of area w^2 / 4k, far below the long double EPS, in all 8 lattice symmetries, flattened by a 4th line
    long long lim = IntHalfplane::LIMIT;
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long k = it % 2 ? lim / 2 - stress::rand_int(0, 1000) : stress::rand_int(1, lim / 2), w = stress::rand_int(1, 2), x0 = stress::rand_int(-lim, lim - 2 * k), y0 = stress::rand_int(-lim, lim - 1);
        int sx = stress::rand_int(0, 1) ? 1 : -1, sy = stress::rand_int(0, 1) ? 1 : -1, swap_xy = stress::rand_int(0, 1);
        auto side = [&](long long ax, long long ay, long long bx, long long by){
            ax *= sx, bx *= sx, ay *= sy, by *= sy;
            if (swap_xy) return sx * sy > 0 ? IntHalfplane(by, bx, ay, ax) : IntHalfplane(ay, ax, by, bx);
            return sx * sy > 0 ? IntHalfplane(ax, ay, bx, by) : IntHalfplane(bx, by, ax, ay);
        };

        vector<IntHalfplane> h = {side(x0, y0, x0 + 1, y0), side(x0 + 2 * k, y0 + 1, x0, y0), side(x0 + w, y0, x0 + w, y0 + 1)};
        if (it % 3 == 0) h.push_back(side(x0 + 1, y0, x0, y0));
        check_exact(h, lim);
    }

    for (long long it = 0; it < stress::scaled(10); it++){
        int n = stress::rand_int(100, 300);
        pair<long long, long long> anchor = {stress::rand_int(-100, 100), stress::rand_int(-100, 100)};
        vector<IntHalfplane> h;
        for (int i = 0; i < n; i++) h.push_back(random_int_halfplane(it % 2 ? 10000 : IntHalfplane::LIMIT, &anchor));
        check_exact(h, IntHalfplane::LIMIT);
    }

    /// The template's own asserts guard the __int128 range and a zero direction
    stress::expect_assert("geometry/halfplane_intersection.cpp", []{ IntHalfplane(0, 0, IntHalfplane::LIMIT + 1, 0); });
    stress::expect_assert("geometry/halfplane_intersection.cpp", []{ IntHalfplane(7, 7, 7, 7); });
    stress::expect_assert("geometry/halfplane_intersection.cpp", []{ halfplane_intersection(vector<IntHalfplane>{}, IntHalfplane::LIMIT + 1); });

    vector<IntHalfplane> big = {IntHalfplane(-1, -1, 1, -1), IntHalfplane(1, -1, 1, 1), IntHalfplane(1, 1, -1, 1), IntHalfplane(-1, 1, -1, -1)};
    while (big.size() < 200000){
        IntHalfplane hp = random_int_halfplane(IntHalfplane::LIMIT, nullptr);
        bool covers = true;
        for (int sx = -1; sx <= 1; sx += 2) for (int sy = -1; sy <= 1; sy += 2) covers &= hp.side(sx, sy) >= 0;
        if (covers) big.push_back(hp);
    }
    shuffle(big.begin(), big.end(), stress::rng());
    vector<Point> square = halfplane_intersection(big, IntHalfplane::LIMIT);
    assert(square.size() == 4 && area2(square) == 8);
}

int main(){
    stress_exact();

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

    /// Lines through one far point: an unbounded cone cut by a huge box, or only that point, which must come out empty
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

    /// Twins of horizontal lines tilted by one ulp across the angle 0 and angle pi seams of the sort
    for (long double y : {-0.3L, 0.3L, -10.0L, 10.0L}){
        long double below = nextafterl(y, y - 1), above = nextafterl(y, y + 1);
        check({Halfplane(Point(0, y), Point(1, below)), Halfplane(Point(5, y), Point(6, y))}, 10);
        check({Halfplane(Point(1, y), Point(0, above)), Halfplane(Point(5, y), Point(4, y))}, 10);
        check({Halfplane(Point(0, y), Point(1, below))}, 10);
        check({Halfplane(Point(1, y), Point(0, above))}, 10);
    }

    /// Antiparallel twins tilted by one ulp around an empty strip: pushing the second next to the first meets parallels
    check({Halfplane(Point(0, nextafterl(-0.6L, -1)), Point(1, -0.6L)),
           Halfplane(Point(1, -0.9L), Point(0, nextafterl(-0.9L, -1))),
           Halfplane(Point(-3, -2), Point(-2, -3))}, 3);

    /// Distinct integer directions with cross 1 must not be merged as one rounded direction
    for (long long k : {30000LL, 1000000LL}){
        check({Halfplane(Point(0, 0), Point(k, k + 1)), Halfplane(Point(0, 0), Point(k - 1, k))}, 1e9);
        check({Halfplane(Point(k, k + 1), Point(0, 0)), Halfplane(Point(k - 1, k), Point(0, 0))}, 1e9);
    }

    /// Decimal 0.1-step inputs repeating lines through other points, so rounded directions are not exactly parallel
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 8), range = it % 2 ? 10 : 100;
        Point anchor = random_point(range);
        vector<Halfplane> h;
        for (int i = 0; i < n; i++){
            Halfplane base = random_halfplane(range, stress::rand_int(0, 1) ? &anchor : nullptr);
            long long dx = llroundl(base.pq.x), dy = llroundl(base.pq.y), g = __gcd(llabs(dx), llabs(dy));
            for (int copies = stress::rand_int(1, 3); copies > 0; copies--){
                int s = stress::rand_int(-3, 3), t = s + stress::rand_int(1, 3);
                Point a = base.p + Point(dx / g, dy / g) * s, b = base.p + Point(dx / g, dy / g) * t;
                h.emplace_back(Point(a.x / 10, a.y / 10), Point(b.x / 10, b.y / 10));
            }
        }
        shuffle(h.begin(), h.end(), stress::rng());
        check(h, it % 3 ? 10 : 1e4);
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
