/***
 *
 * Pick's Thoerem and other related methods
 * Calculates area of polygon and number of lattice points inside and on border
 *
 * Complexity:
 *   - area2: O(n)
 *   - on_border, on_interior: O(n log C), one gcd per edge, C the largest coordinate difference along an edge
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

/// number of lattice points strictly on the polygon border (edges)
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
    auto polygon = {
        Point(0, 0), Point(10, 0), Point(10, 10), Point(5, 20), Point(0, 10)
    };

    assert(area2(polygon) == 300);
    assert(on_border(polygon) == 40);
    assert(on_interior(polygon) == 131);

    return 0;
}
