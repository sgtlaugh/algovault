#include "../common.h"

#define main library_main
#include "../../code_library/geometry/polygon_union.cpp"
#undef main

typedef vector<vector<Point>> Polygons;

long long turn(Point o, Point a, Point b){
    return (long long)((a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x));
}

bool on_segment(Point p, Point a, Point b){
    return turn(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

bool segments_touch(Point a, Point b, Point c, Point d){
    auto sgn = [](long long v){ return (v > 0) - (v < 0); };
    if (sgn(turn(a, b, c)) * sgn(turn(a, b, d)) < 0 && sgn(turn(c, d, a)) * sgn(turn(c, d, b)) < 0) return true;
    return on_segment(c, a, b) || on_segment(d, a, b) || on_segment(a, c, d) || on_segment(b, c, d);
}

/// Simple: adjacent edges share only their common vertex, non-adjacent edges are disjoint, area is not zero
bool is_simple(const vector<Point>& poly){
    int n = poly.size();
    long long area = 0;
    for (int i = 0; i < n; i++) area += turn(poly[0], poly[i], poly[(i + 1) % n]);
    if (n < 3 || area == 0) return false;

    for (int i = 0; i < n; i++){
        Point a = poly[i], b = poly[(i + 1) % n], c = poly[(i + 2) % n];
        if (on_segment(a, b, c) || on_segment(c, a, b)) return false;
        for (int j = i + 2; j < n; j++){
            if (i == 0 && j == n - 1) continue;
            if (segments_touch(a, b, poly[j], poly[(j + 1) % n])) return false;
        }
    }
    return true;
}

Point random_point(int lo, int hi){
    return {(Float)stress::rand_int(lo, hi), (Float)stress::rand_int(lo, hi)};
}

/// Star-shaped around a random center, retried until simple: covers convex, concave and collinear-vertex polygons
vector<Point> random_polygon(int range){
    while (true){
        int kind = stress::rand_int(0, 3), k = stress::rand_int(3, 7);
        vector<Point> poly;
        if (kind == 0){
            Point p = random_point(0, range - 1), q = random_point(0, range - 1);
            Float x1 = min(p.x, q.x), y1 = min(p.y, q.y), x2 = max(p.x, q.x) + 1, y2 = max(p.y, q.y) + 1;
            poly = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
        }
        else{
            if (kind == 1) k = 3;
            Point center = random_point(0, range);
            for (int i = 0; i < k; i++) poly.push_back(random_point(0, range));
            sort(poly.begin(), poly.end(), [&](const Point& p, const Point& q){
                return atan2l(p.y - center.y, p.x - center.x) < atan2l(q.y - center.y, q.x - center.x);
            });
        }

        if (!is_simple(poly)) continue;
        if (stress::rand_int(0, 1)) reverse(poly.begin(), poly.end());
        rotate(poly.begin(), poly.begin() + stress::rand_int(0, poly.size() - 1), poly.end());
        return poly;
    }
}

/// Zero-area input the header says adds 0: a segment, or three distinct collinear points in any order
vector<Point> random_degenerate(int range){
    if (stress::rand_int(0, 1)) return {random_point(0, range), random_point(0, range)};
    while (true){
        Point p = random_point(0, range), d = random_point(-2, 2);
        if (d.x == 0 && d.y == 0) continue;
        vector<Point> poly = {p, {p.x + d.x, p.y + d.y}, {p.x + 2 * d.x, p.y + 2 * d.y}};
        if (poly[2].x < 0 || poly[2].x > range || poly[2].y < 0 || poly[2].y > range) continue;
        shuffle(poly.begin(), poly.end(), stress::rng());
        return poly;
    }
}

Polygons random_polygons(int range, int count){
    Polygons polys;
    for (int i = 0; i < count; i++){
        if (!polys.empty() && stress::rand_int(0, 5) == 0){
            vector<Point> copy = polys[stress::rand_int(0, polys.size() - 1)];
            if (stress::rand_int(0, 1)) reverse(copy.begin(), copy.end());
            rotate(copy.begin(), copy.begin() + stress::rand_int(0, copy.size() - 1), copy.end());
            polys.push_back(copy);
        }
        else polys.push_back(random_polygon(range));
    }
    return polys;
}

/// Vertical strips cut at every vertex and edge crossing: inside a strip no edges cross, the covered length is linear in x and its midpoint value is exact
Float strip_union(const Polygons& polys){
    vector<pair<Point, Point>> edges;
    for (auto& poly: polys){
        for (int i = 0; i < (int)poly.size(); i++) edges.push_back({poly[i], poly[(i + 1) % poly.size()]});
    }

    vector<Float> xs;
    for (auto& [p, q]: edges) xs.push_back(p.x);
    for (auto& [p, q]: edges){
        for (auto& [r, s]: edges){
            Float rx = q.x - p.x, ry = q.y - p.y, sx = s.x - r.x, sy = s.y - r.y, den = rx * sy - ry * sx;
            if (den == 0) continue;
            Float t = ((r.x - p.x) * sy - (r.y - p.y) * sx) / den, w = ((r.x - p.x) * ry - (r.y - p.y) * rx) / den;
            if (t >= 0 && t <= 1 && w >= 0 && w <= 1) xs.push_back(p.x + t * rx);
        }
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());

    Float res = 0;
    for (int s = 0; s + 1 < (int)xs.size(); s++){
        Float mid = (xs[s] + xs[s + 1]) / 2;
        vector<pair<Float, Float>> intervals;
        for (auto& poly: polys){
            vector<Float> ys;
            for (int i = 0; i < (int)poly.size(); i++){
                Point p = poly[i], q = poly[(i + 1) % poly.size()];
                if ((p.x < mid) != (q.x < mid)) ys.push_back(p.y + (q.y - p.y) * (mid - p.x) / (q.x - p.x));
            }
            sort(ys.begin(), ys.end());
            for (int i = 0; i + 1 < (int)ys.size(); i += 2) intervals.push_back({ys[i], ys[i + 1]});
        }

        sort(intervals.begin(), intervals.end());
        Float covered = 0, reach = -1e30;
        for (auto& [lo, hi]: intervals){
            if (hi <= reach) continue;
            covered += hi - max(lo, reach);
            reach = hi;
        }
        res += covered * (xs[s + 1] - xs[s]);
    }
    return res;
}

/// Even-odd point test on a sample grid, error bounded by the sample cells each edge passes through
void check_grid(const Polygons& polys, int range, Float expected){
    const int cells = 300;
    Float h = (Float)range / cells, bound = 0;
    for (auto& poly: polys){
        for (int i = 0; i < (int)poly.size(); i++){
            Point d = poly[(i + 1) % poly.size()] - poly[i];
            bound += (fabsl(d.x) + fabsl(d.y) + 2 * h) * h;
        }
    }

    long long inside = 0;
    for (int gx = 0; gx < cells; gx++){
        for (int gy = 0; gy < cells; gy++){
            Float x = (gx + 0.5L) * h, y = (gy + 0.5L) * h;
            bool any = false;
            for (auto& poly: polys){
                bool in = false;
                for (int i = 0, n = poly.size(); i < n && !any; i++){
                    Point p = poly[i], q = poly[(i + 1) % n];
                    if ((p.y > y) != (q.y > y) && x < p.x + (q.x - p.x) * (y - p.y) / (q.y - p.y)) in = !in;
                }
                any = any || in;
            }
            inside += any;
        }
    }
    assert(fabsl(inside * h * h - expected) <= bound);
}

Polygons transform(const Polygons& polys, Float scale, Float dx, Float dy){
    Polygons res = polys;
    for (auto& poly: res){
        for (auto& p: poly) p = {p.x * scale + dx, p.y * scale + dy};
    }
    return res;
}

/// Integer rectangles on a 100 x 100 grid, union counted cell by cell
void check_rectangles(int count){
    Polygons polys;
    vector<vector<bool>> covered(100, vector<bool>(100, false));
    for (int i = 0; i < count; i++){
        int x1 = stress::rand_int(0, 99), y1 = stress::rand_int(0, 99);
        int x2 = stress::rand_int(x1 + 1, min(100, x1 + 20)), y2 = stress::rand_int(y1 + 1, min(100, y1 + 20));
        polys.push_back({{(Float)x1, (Float)y1}, {(Float)x2, (Float)y1}, {(Float)x2, (Float)y2}, {(Float)x1, (Float)y2}});
        for (int x = x1; x < x2; x++){
            for (int y = y1; y < y2; y++) covered[x][y] = true;
        }
    }

    long long cells = 0;
    for (auto& row: covered) cells += count_if(row.begin(), row.end(), [](bool c){ return c; });
    assert(fabsl(polygon_union(polys) - cells) <= 1e-9L * max(1LL, cells));
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int range = it % 3 == 0 ? 3 : (it % 3 == 1 ? 8 : 30), count = stress::rand_int(1, it % 5 == 0 ? 8 : 4);
        Polygons polys = random_polygons(range, count);
        if (stress::rand_int(0, 4) == 0) polys.insert(polys.begin() + stress::rand_int(0, polys.size()), random_degenerate(range));

        Float expected = strip_union(polys), res = polygon_union(polys);
        assert(fabsl(res - expected) <= 1e-9L * max((Float)1, expected));

        Polygons shuffled = polys;
        shuffle(shuffled.begin(), shuffled.end(), stress::rng());
        assert(fabsl(polygon_union(shuffled) - expected) <= 1e-9L * max((Float)1, expected));

        /// Header boundary: integer coordinates reaching |x|, |y| = 1e9, translated to the corner and scaled across the whole range
        Polygons cornered = transform(polys, 1, 1e9 - range, -1e9);
        assert(fabsl(polygon_union(cornered) - expected) <= 1e-6L * max((Float)1, expected));
        Float scale = floorl(2e9L / range);
        Polygons spread = transform(polys, scale, -1e9, -1e9);
        assert(fabsl(polygon_union(spread) / (scale * scale) - expected) <= 1e-12L * max((Float)1, expected));

        if (it % 100 == 0 && range <= 8) check_grid(polys, range, expected);
    }

    for (long long it = 0; it < stress::scaled(20); it++) check_rectangles(stress::rand_int(1, 60));
    check_rectangles(500);

    return 0;
}
