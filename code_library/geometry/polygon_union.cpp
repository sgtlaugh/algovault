/***
 *
 * Polygon Union Area
 * Area of the union of n simple polygons, convex or not, in any orientation
 *
 * Complexity: O(N^2 log N), N = total number of vertices over all polygons
 *
 * polygon_union(polys): each polygon is a list of distinct vertices in boundary order, clockwise or counter-clockwise
 *     polygons may overlap, nest, touch and share edges; polygons with fewer than 3 vertices or all vertices collinear add 0
 *
 * For every edge, the parts covered by other polygons are cut out and the area is the shoelace sum over the visible parts
 * Collinear edges going the same way are covered only by the polygon with the lower index, so duplicates count once
 *
 * Sign tests are exact for integer coordinates with |x|, |y| <= 1e9 in 80-bit long double (g++ on x86)
 * Fractional coordinates work, but a shared edge may then be missed or double counted by rounding
 *
 * Example:
 *   polygon_union({{{0, 0}, {2, 0}, {2, 2}, {0, 2}}, {{1, 1}, {3, 1}, {3, 3}, {1, 3}}})  // 7
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

typedef long double Float;

struct Point{
    Float x, y;

    Point operator-(const Point& p) const{
        return {x - p.x, y - p.y};
    }

    Float cross(const Point& p) const{
        return x * p.y - y * p.x;
    }

    Float dot(const Point& p) const{
        return x * p.x + y * p.y;
    }
};

int sign(Float v){
    return (v > 0) - (v < 0);
}

/// Turn at the lowest-leftmost vertex: both neighbours lie on one side of it, so the sign is exact and 0 only for degenerate polygons
int orientation(const vector<Point>& poly){
    int n = poly.size(), k = 0;
    if (n < 3) return 0;

    for (int i = 1; i < n; i++){
        if (poly[i].x < poly[k].x || (poly[i].x == poly[k].x && poly[i].y < poly[k].y)) k = i;
    }
    return sign((poly[(k + 1) % n] - poly[k]).cross(poly[(k + n - 1) % n] - poly[k]));
}

Float polygon_union(vector<vector<Point>> polys){
    for (auto& poly: polys){
        if (orientation(poly) < 0) reverse(poly.begin(), poly.end());
    }
    polys.erase(remove_if(polys.begin(), polys.end(), [](const vector<Point>& poly){ return orientation(poly) == 0; }), polys.end());

    auto ratio = [](const Point& p, const Point& d){ return sign(d.x) ? p.x / d.x : p.y / d.y; };
    Float res = 0;
    int m = polys.size();
    for (int i = 0; i < m; i++){
        int n = polys[i].size();
        for (int v = 0; v < n; v++){
            Point a = polys[i][v], b = polys[i][(v + 1) % n], ab = b - a;
            vector<pair<Float, int>> events = {{0, 0}, {1, 0}};

            for (int j = 0; j < m; j++){
                if (i == j) continue;
                int k = polys[j].size();
                for (int u = 0; u < k; u++){
                    Point c = polys[j][u], d = polys[j][(u + 1) % k];
                    int sc = sign(ab.cross(c - a)), sd = sign(ab.cross(d - a));
                    if (sc != sd){
                        Float sa = (d - c).cross(a - c), sb = (d - c).cross(b - c);
                        if (min(sc, sd) < 0) events.push_back({sa / (sa - sb), sign(sc - sd)});
                    }
                    else if (sc == 0 && j < i && sign(ab.dot(d - c)) > 0){
                        events.push_back({ratio(c - a, ab), 1});
                        events.push_back({ratio(d - a, ab), -1});
                    }
                }
            }

            sort(events.begin(), events.end());
            for (auto& e: events) e.first = min(max(e.first, (Float)0), (Float)1);

            Float visible = 0;
            int cover = events[0].second;
            for (int e = 1; e < (int)events.size(); e++){
                if (!cover) visible += events[e].first - events[e - 1].first;
                cover += events[e].second;
            }
            res += a.cross(b) * visible;
        }
    }
    return res / 2;
}

int main(){
    vector<Point> square = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    vector<Point> shifted = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
    vector<Point> clockwise = {{0, 0}, {0, 2}, {2, 2}, {2, 0}};
    vector<Point> rotated = {{2, 2}, {0, 2}, {0, 0}, {2, 0}};

    assert(abs(polygon_union({}) - 0) < 1e-9);
    assert(abs(polygon_union({square}) - 4) < 1e-9);
    assert(abs(polygon_union({clockwise}) - 4) < 1e-9);
    assert(abs(polygon_union({square, shifted}) - 7) < 1e-9);
    assert(abs(polygon_union({square, square}) - 4) < 1e-9);
    assert(abs(polygon_union({square, rotated, clockwise}) - 4) < 1e-9);
    assert(abs(polygon_union({{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, {{1, 0}, {2, 0}, {2, 1}, {1, 1}}}) - 2) < 1e-9);
    assert(abs(polygon_union({{{0, 0}, {4, 0}, {4, 4}, {0, 4}}, {{1, 1}, {3, 1}, {1, 3}}}) - 16) < 1e-9);
    assert(abs(polygon_union({{{0, 0}, {2, 0}, {0, 2}}, {{5, 5}, {7, 5}, {5, 7}}}) - 4) < 1e-9);
    assert(abs(polygon_union({square, {{1, 1}, {3, 1}, {1, 3}}}) - 5) < 1e-9);

    vector<Point> l_shape = {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}};
    assert(abs(polygon_union({l_shape}) - 3) < 1e-9);
    assert(abs(polygon_union({l_shape, {{1, 1}, {2, 1}, {2, 2}, {1, 2}}}) - 4) < 1e-9);
    vector<Point> raised = {{0, 1}, {2, 1}, {2, 3}, {0, 3}};
    assert(abs(polygon_union({{{0, 1}, {5, 1}}, {{2, 0}, {2, 2}, {2, 4}}, raised}) - 4) < 1e-9);

    const Float far = 1e9;
    vector<Point> corner = {{far - 2, -far}, {far, -far}, {far, -far + 2}, {far - 2, -far + 2}};
    vector<Point> near_corner = {{far - 3, -far + 1}, {far - 1, -far + 1}, {far - 1, -far + 3}, {far - 3, -far + 3}};
    assert(abs(polygon_union({corner, near_corner}) - 7) < 1e-9);
    assert(abs(polygon_union({{{-far, -far}, {far, -far}, {far, far}, {-far, far}}, corner}) - 4e18) < 1e3);

    return 0;
}
