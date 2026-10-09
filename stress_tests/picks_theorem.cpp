#include "common.h"

#define main library_main
#include "../code_library/picks_theorem.cpp"
#undef main

/// Vertices sorted by angle around a centre that every angular gap keeps strictly inside, so the polygon is simple
vector<Point> star_polygon(int range){
    while (true){
        long long cx = stress::rand_int(-range, range), cy = stress::rand_int(-range, range);
        map<long double, Point> by_angle;
        for (int k = stress::rand_int(3, 12); k; k--){
            Point p(stress::rand_int(-range, range), stress::rand_int(-range, range));
            if (p.x == cx && p.y == cy) continue;
            by_angle.emplace(atan2l(p.y - cy, p.x - cx), p);
        }
        if (by_angle.size() < 3) continue;

        vector<pair<long double, Point>> v(by_angle.begin(), by_angle.end());
        bool inside = true;
        for (size_t i = 0; i < v.size(); i++){
            long double next = i + 1 < v.size() ? v[i + 1].first : v[0].first + 2 * acosl(-1);
            Point a(v[i].second.x - cx, v[i].second.y - cy), b(v[(i + 1) % v.size()].second.x - cx, v[(i + 1) % v.size()].second.y - cy);
            if (next - v[i].first >= acosl(-1) || a.x * b.y - a.y * b.x <= 0) inside = false;
        }
        if (!inside) continue;

        vector<Point> poly;
        for (auto& [angle, p] : v) poly.push_back(p);
        if (stress::rand_int(0, 1)) reverse(poly.begin(), poly.end());
        rotate(poly.begin(), poly.begin() + stress::rand_int(0, poly.size() - 1), poly.end());
        return poly;
    }
}

int main(){
    /// A sliver far from the y-axis, its shoelace terms reach 3.2e19 though twice its area is only 4e9
    vector<Point> sliver = {Point(4000000000LL, -2000000000LL), Point(4000000001LL, -2000000000LL), Point(4000000000LL, 2000000000LL)};
    assert(area2(sliver) == 4000000000LL && on_border(sliver) == 4000000002LL && on_interior(sliver) == 0);

    for (long long it = 0; it < stress::scaled(4000); it++){
        int range = stress::rand_int(2, it % 5 ? 8 : 40);
        auto poly = star_polygon(range);
        int n = poly.size();

        long long border = 0, interior = 0;
        for (long long x = -range; x <= range; x++){
            for (long long y = -range; y <= range; y++){
                bool on_edge = false, odd = false;
                for (int i = 0, j = n - 1; i < n; j = i++){
                    const Point &a = poly[j], &b = poly[i];
                    if ((b.x - a.x) * (y - a.y) == (b.y - a.y) * (x - a.x) && min(a.x, b.x) <= x && x <= max(a.x, b.x) && min(a.y, b.y) <= y && y <= max(a.y, b.y)) on_edge = true;
                    if ((a.y > y) != (b.y > y)){
                        /// Exact crossing test of the rightward ray, x_cross > x without division
                        long long num = (b.x - a.x) * (y - a.y) + a.x * (b.y - a.y) - x * (b.y - a.y);
                        if ((b.y - a.y > 0) ? num > 0 : num < 0) odd = !odd;
                    }
                }
                if (on_edge) border++;
                else if (odd) interior++;
            }
        }

        assert(on_border(poly) == border);
        assert(on_interior(poly) == interior);
        assert(area2(poly) == 2 * interior + border - 2);
    }

    /// Large coordinates: right triangles with legs along the axes, where every count has a closed form
    for (long long it = 0; it < stress::scaled(1000); it++){
        long long x0 = stress::rand_int(-1000000000, 1000000000), y0 = stress::rand_int(-1000000000, 1000000000);
        long long w = stress::rand_int(1, 1000000000), h = stress::rand_int(1, 1000000000), g = __gcd(w, h);
        vector<Point> tri = {Point(x0, y0), Point(x0 + w, y0), Point(x0, y0 + h)};
        if (stress::rand_int(0, 1)) reverse(tri.begin(), tri.end());

        long long border = w + h + g;
        assert(area2(tri) == w * h);
        assert(on_border(tri) == border);
        assert(on_interior(tri) == (w * h - border + 2) / 2);
    }

    /// Large star polygons, checked against a triangle fan in __int128
    for (long long it = 0; it < stress::scaled(1000); it++){
        auto poly = star_polygon(1000000000);
        __int128 fan = 0;
        for (size_t i = 1; i + 1 < poly.size(); i++){
            __int128 ax = poly[i].x - poly[0].x, ay = poly[i].y - poly[0].y, bx = poly[i + 1].x - poly[0].x, by = poly[i + 1].y - poly[0].y;
            fan += ax * by - ay * bx;
        }
        assert(area2(poly) == (long long)(fan < 0 ? -fan : fan));
    }
    return 0;
}
