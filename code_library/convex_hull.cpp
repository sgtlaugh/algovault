/***
 *
 * Convex hull using the Monotone Chain algorithm
 *
 * Complexity: O(N log N)
 *
 * Optimization notes: can be converted to O(n) using radix sort
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Point {
    int64_t x, y; /// x*x or y*y should not overflow

    Point(){}
    Point(int64_t x, int64_t y) : x(x), y(y) {}

    inline bool operator < (const Point &p) const {
        return ((x < p.x) || (x == p.x && y < p.y));
    }
};

int64_t cross(const Point &O, const Point &A, const Point &B){
    return ((A.x - O.x) * (B.y - O.y)) - ((A.y - O.y) * (B.x - O.x));
}

/***
 *
 * Returns the strict convex hull in counter-clockwise order, starting from the lowest x (then lowest y) point
 * Collinear and duplicate points are dropped, all-collinear input gives the two endpoints
 *
***/

vector<Point> get_convex_hull(vector<Point> P){
    sort(P.begin(), P.end());
    P.erase(unique(P.begin(), P.end(), [](const Point& a, const Point& b){ return a.x == b.x && a.y == b.y; }), P.end());

    int i, t, k = 0, n = P.size();
    if (n <= 1) return P;
    vector<Point> H(n << 1);

    for (i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }
    for (i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    H.resize(k - 1);
    return H;
}

/***
 *
 * Returns whether the polygon is convex or not
 * Points in P are given in clockwise or anti-clockwise order, collinear vertices are allowed
 *
***/

bool is_convex(const vector <Point>& P){
    int n = P.size(), sign = 0, flips = 0, first_dy = 0, last_dy = 0;
    if (n <= 2) return false; /// Line or point is not convex

    for (int i = 0; i < n; i++){
        const Point &a = P[i], &b = P[(i + 1) % n], &c = P[(i + 2) % n];
        int64_t turn = cross(a, b, c);
        if (turn){
            if (sign && (turn > 0) != (sign > 0)) return false;
            sign = turn > 0 ? 1 : -1;
        }

        int dy = (b.y > a.y) - (b.y < a.y);
        if (dy){
            if (last_dy && dy != last_dy) flips++;
            if (!first_dy) first_dy = dy;
            last_dy = dy;
        }
    }
    if (first_dy != last_dy) flips++;

    /// Turning one way everywhere still allows star polygons, winding exactly once means y changes direction twice
    return sign != 0 && flips == 2;
}

int main(){
    vector <Point> polygon = {Point(0, 0), Point(0, 10), Point(1, 1), Point(2, 20), Point(5, 5), Point(10, 10), Point(10, 0)};
    assert(!is_convex(polygon));

    vector <Point> hull = get_convex_hull(polygon);
    assert(is_convex(hull));

    vector <Point> expected_hull = {Point(0, 0), Point(10, 0), Point(10, 10), Point(2, 20), Point(0, 10)};
    assert((int)hull.size() == (int)expected_hull.size());

    for (int i = 0; i < (int)hull.size(); i++){
        assert(hull[i].x == expected_hull[i].x && hull[i].y == expected_hull[i].y);
    }

    return 0;
}
