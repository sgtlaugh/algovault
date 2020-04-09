#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Point {
    long long x, y; /// x*x or y*y should fit long long because of cross() function

    Point(){}
    Point(long long x, long long y) : x(x), y(y) {}

    inline bool operator < (const Point &p) const {
        return ((x < p.x) || (x == p.x && y < p.y));
    }
};

long long cross(const Point &O, const Point &A, const Point &B){
    return ((A.x - O.x) * (B.y - O.y)) - ((A.y - O.y) * (B.x - O.x));
}

/// Returns the convex hull using the monotone chain algorithm
vector <Point> get_convex_hull(vector<Point> P){
    int i, t, k = 0, n = P.size();
    vector<Point> H(n << 1);
    sort(P.begin(), P.end()); /// Can be converted to O(n) using radix sort

    for (i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) < 0) k--;
        H[k++] = P[i];
    }
    for (i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(H[k - 2], H[k - 1], P[i]) < 0) k--;
        H[k++] = P[i];
    }

    H.resize(k - 1); /// The last point is the same as the first in this implementation
    return H;
}

/***
 *
 * returns whether the polygon is convex or not
 * points in P is given in clock-wise of anti-clockwise order
 *
***/

bool is_convex(vector <Point> P){
    int n = P.size();
    if (n <= 2) return false; /// Line or point is not convex

    n++, P.push_back(P[0]);  /// Last point = first point
    bool flag1 = (cross(P[0], P[1], P[2]) > 0);
    for (int i = 1; (i + 1) < n; i++){
        bool flag2 = (cross(P[i], P[i + 1], (i + 2) == n ? P[1] : P[i + 2]) > 0);
        if (flag1 ^ flag2) return false;
    }

    return true;
}

int main(){
    vector <Point> polygon = {Point(0, 0), Point(0, 10), Point(1, 1), Point(2, 20), Point(5, 5), Point(10, 10), Point(10, 0)};
    printf("%d\n", is_convex(polygon));               /// False

    vector <Point> hull = get_convex_hull(polygon);
    printf("%d\n", is_convex(hull));                  /// True

    printf("Convex hull contains %d points\n", (int)hull.size());
    for (auto point: hull){
        printf("(%lld, %lld), ", point.x, point.y);   /// (0, 0), (10, 0), (10, 10), (2, 20), (0, 10)
    }
    puts("");

    return 0;
}
