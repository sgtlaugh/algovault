/***
 *
 * Pick's Theorem and other related methods
 * Calculates area of polygon and number of lattice points inside and on border
 *
 * Complexity:
 *   - area2: O(n)
 *   - on_border, on_interior: O(n log C), one gcd per edge, C the largest coordinate difference along an edge
 *
 * poly is a simple polygon with lattice vertices in either orientation, |x|, |y| <= 1e9
 *   area2(poly): twice the area, exact
 *   on_border(poly): the count of lattice points on the border, vertices included
 *   on_interior(poly): the count of lattice points strictly inside, from Pick's theorem A = I + B / 2 - 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;

    Point() {}
    Point(long long x, long long y) : x(x), y(y) {}
};

/// area of the polygon multiplied by 2
/// Partial sums overflow long long near 1e9 coordinates even when the area fits
long long area2(const vector<Point>& poly){
    __int128 res = 0;
    int i, j, n = poly.size();
    for (i = 0, j = n - 1; i < n; j = i++){
        res += (__int128)(poly[j].x + poly[i].x) * (poly[j].y - poly[i].y);
    }
    return res < 0 ? -res : res;
}

/// number of lattice points on the polygon border, vertices included
long long on_border(const vector<Point>& poly){
    long long res = 0;
    int i, j, n = poly.size();
    for (i = 0, j = n - 1; i < n; j = i++){
        res += __gcd(abs(poly[i].x - poly[j].x), abs(poly[i].y - poly[j].y));
    }
    return res;
}

/// number of lattice points strictly inside the polygon
long long on_interior(const vector<Point>& poly){
    long long res = 2 + area2(poly) - on_border(poly);
    return res >> 1;
}

int main(){
    int n, m, cs = 0;

    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        vector<Point> pts(n), poly;
        for (auto& p: pts) scanf("%lld %lld", &p.x, &p.y);

        for (int i = 0; i < n; i++){
            if (area2({pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]})) poly.push_back(pts[i]);
        }

        long long g = 0;
        Point o = poly[0];
        for (auto& p: poly) g = __gcd(g, __gcd(abs(p.x - o.x), abs(p.y - o.y)));
        for (auto& p: poly) p = Point((p.x - o.x) / g, (p.y - o.y) / g);

        long long a = area2(poly), b = on_border(poly), res = on_interior(poly);
        for (long long i = 2; i <= m; i++) res += (a * i * i - b * i + 2) / 2;
        printf("Case %d: %lld\n", ++cs, res);
    }
    return 0;
}
