/***
 *
 * Closest Pair of Points
 * Finds two points with the minimum Euclidean distance, using exact integer squared distances
 *
 * Complexity: O(n log n) time, O(n) space
 *
 * closest_pair(points) returns {dist2, i, j}: the squared distance and the indices i < j of one closest pair
 * Requires n >= 2 and |x|, |y| <= 1e9, so a squared distance (at most 8e18) fits in a long long
 * Duplicate points are allowed and give dist2 = 0
 *
 * Sweeps the points by y, keeping the ones within the current best distance in a set ordered by x
 * Each new point only checks the set entries inside a d x 2d box, which holds O(1) points
 *
 * Example:
 *   ClosestPairResult res = closest_pair({{0, 0}, {5, 5}, {1, 2}});
 *   // res.dist2 = 5, res.i = 0, res.j = 2
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;
};

struct ClosestPairResult{
    long long dist2;
    int i, j;
};

ClosestPairResult closest_pair(const vector<Point>& points){
    int n = points.size();
    assert(n >= 2);

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){ return points[a].y < points[b].y; });

    ClosestPairResult res = {LLONG_MAX, -1, -1};
    set<tuple<long long, long long, int>> window;
    for (int k = 0, tail = 0; k < n && res.dist2 > 0; k++){
        const Point& p = points[order[k]];
        /// Strictly above the true distance even when sqrtl rounds, so no closer point is pruned
        long long d = (long long)sqrtl((long double)res.dist2) + 1;
        while (points[order[tail]].y <= p.y - d){
            const Point& q = points[order[tail]];
            window.erase({q.x, q.y, order[tail++]});
        }

        auto it = window.lower_bound({p.x - d + 1, LLONG_MIN, INT_MIN});
        for (; it != window.end() && get<0>(*it) < p.x + d; it++){
            long long dx = get<0>(*it) - p.x, dy = get<1>(*it) - p.y;
            if (dx * dx + dy * dy < res.dist2) res = {dx * dx + dy * dy, get<2>(*it), order[k]};
        }
        window.insert({p.x, p.y, order[k]});
    }

    if (res.i > res.j) swap(res.i, res.j);
    return res;
}

int main(){
    auto check = [](const vector<Point>& points, long long dist2, int i, int j){
        ClosestPairResult res = closest_pair(points);
        return res.dist2 == dist2 && res.i == i && res.j == j;
    };

    assert(check({{0, 0}, {3, 4}}, 25, 0, 1));
    assert(check({{0, 0}, {5, 5}, {1, 2}}, 5, 0, 2));
    assert(check({{2, 3}, {12, 30}, {40, 50}, {5, 1}, {12, 10}, {3, 4}}, 2, 0, 5));
    assert(check({{0, 0}, {0, 7}, {0, 3}, {0, 12}}, 9, 0, 2));
    assert(check({{-4, 9}, {6, 9}, {21, 9}, {13, 9}}, 49, 1, 3));
    assert(check({{5, 5}, {1, 2}, {5, 5}}, 0, 0, 2));
    assert(check({{-1000000000, -1000000000}, {1000000000, 1000000000}}, 8000000000000000000LL, 0, 1));
    assert(check({{-1000000000, 1000000000}, {1000000000, -1000000000}, {1000000000, 999999999}}, 3999999996000000001LL, 1, 2));

    return 0;
}
