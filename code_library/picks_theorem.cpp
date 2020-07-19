/***
 *
 * Pick's Thoerem and other related methods
 * Calculates area of polygon in O(n)
 * Calculates number of lattice points inside and on border in O(n log max(P_i))
 *
*/

#include <bits/stdc++.h>

using namespace std;

struct Point{
    long long x, y;
    Point() {}
    Point(long long x, long long y) : x(x), y(y) {}
};

/// area of the polygon multiplied by 2
long long area2(const vector<Point>& poly){
    long long res = 0;
    int i, j, n = poly.size();
    for (i = 0, j = n - 1; i < n; j = i++){
        res += ((poly[j].x + poly[i].x) * (poly[j].y - poly[i].y));
    }
    return abs(res);
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

    printf("%lld\n", area2(polygon));        /// 300
    printf("%lld\n", on_border(polygon));    /// 40
    printf("%lld\n", on_interior(polygon));  /// 131

    return 0;
}
