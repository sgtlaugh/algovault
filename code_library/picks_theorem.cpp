#include <bits/stdc++.h>

using namespace std;

struct Point{
    long long x, y;
    Point() {}
    Point(long long x, long long y) : x(x), y(y) {}
};

/// returns the area of the polygon multiplied by 2
long long area2(int n, Point poly[]){
    long long res = 0;
    for (int i = 0, j = n - 1; i < n; j = i++){
        res += ((poly[j].x + poly[i].x) * (poly[j].y - poly[i].y));
    }
    return abs(res);
}

/// number of lattice points strictly on the polygon border (edges)
long long on_border(int n, Point poly[]){
    long long res = 0;
    for (int i = 0, j = n - 1; i < n; j = i++){
        res += __gcd(abs(poly[i].x - poly[j].x), abs(poly[i].y - poly[j].y));
    }
    return res;
}

/// number of lattice points strictly inside polygon
long long on_interior(int n, Point poly[]){
    long long res = 2 + area2(n, poly) - on_border(n, poly);
    return res >> 1;
}

int main(){
    Point polygon[] = {
        Point(0, 0), Point(10, 0), Point(10, 10), Point(5, 20), Point(0, 10)
    };

    printf("%lld\n", area2(5, polygon));        /// 300
    printf("%lld\n", on_border(5, polygon));    /// 40
    printf("%lld\n", on_interior(5, polygon));  /// 131

    return 0;
}
