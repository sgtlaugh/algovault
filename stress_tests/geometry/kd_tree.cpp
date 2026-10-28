#include "../common.h"

#define main library_main
#include "../../code_library/geometry/kd_tree.cpp"
#undef main

const long long C = 1000000000;

long long brute_dist(const Point& a, const Point& b, KDTree::Metric metric){
    long long dx = a.x - b.x, dy = a.y - b.y;
    return metric == KDTree::EUCLIDEAN ? dx * dx + dy * dy : abs(dx) + abs(dy);
}

void check_nearest(const vector<Point>& points, const KDTree& tree, const Point& q, int skip){
    long long best = LLONG_MAX;
    for (int j = 0; j < (int)points.size(); j++){
        if (j != skip) best = min(best, brute_dist(points[j], q, tree.metric));
    }

    KDTree::Neighbor res = tree.nearest(q, skip);
    assert(res.dist == best);
    if (best == LLONG_MAX){
        assert(res.i == -1);
        return;
    }
    assert(0 <= res.i && res.i < (int)points.size() && res.i != skip);
    assert(brute_dist(points[res.i], q, tree.metric) == best);
}

void check_report(const vector<Point>& points, const KDTree& tree, long long x1, long long y1, long long x2,
                  long long y2){
    vector<int> expected;
    for (int j = 0; j < (int)points.size(); j++){
        if (x1 <= points[j].x && points[j].x <= x2 && y1 <= points[j].y && points[j].y <= y2) expected.push_back(j);
    }

    vector<int> res = tree.report(x1, y1, x2, y2);
    sort(res.begin(), res.end());
    assert(res == expected);
}

/// Often one value wide (catches a strict < on a box edge), sometimes stretched to the coordinate bound
pair<long long, long long> random_interval(long long range){
    long long lo = stress::rand_int(-range, range), hi = stress::rand_int(-range, range);
    if (stress::rand_int(0, 4) == 0) hi = lo;
    if (stress::rand_int(0, 9) == 0) lo = -C;
    if (stress::rand_int(0, 9) == 0) hi = C;

    return {min(lo, hi), max(lo, hi)};
}

/// Every input point with itself skipped (the FAILURE query) plus random queries and rectangles around the set
void check_all(const vector<Point>& points, int queries, long long range){
    KDTree euclidean(points);
    for (int k = 0; k < queries; k++){
        auto [x1, x2] = random_interval(range);
        auto [y1, y2] = random_interval(range);
        check_report(points, euclidean, x1, y1, x2, y2);
    }

    for (KDTree::Metric metric : {KDTree::EUCLIDEAN, KDTree::MANHATTAN}){
        KDTree tree(points, metric);
        for (int i = 0; i < (int)points.size(); i++) check_nearest(points, tree, points[i], i);
        for (int k = 0; k < queries; k++){
            Point q = {stress::rand_int(-range, range), stress::rand_int(-range, range)};
            check_nearest(points, tree, q, stress::rand_int(-1, (int)points.size() - 1));
        }
    }
}

/// All multisets of up to 3 points on a 3 x 2 grid, every query point and rectangle near the grid and every skip
void check_exhaustive(){
    vector<Point> cells;
    for (int x = 0; x < 3; x++){
        for (int y = 0; y < 2; y++) cells.push_back({x, y});
    }

    for (int n = 0; n <= 3; n++){
        int total = 1;
        for (int k = 0; k < n; k++) total *= cells.size();
        for (int mask = 0; mask < total; mask++){
            vector<Point> points;
            int rest = mask;
            for (int k = 0; k < n; k++){
                points.push_back(cells[rest % cells.size()]);
                rest /= cells.size();
            }
            for (KDTree::Metric metric : {KDTree::EUCLIDEAN, KDTree::MANHATTAN}){
                KDTree tree(points, metric);
                for (int x1 = -1; x1 <= 3; x1++){
                    for (int x2 = x1; x2 <= 3; x2++){
                        for (int y1 = -1; y1 <= 2; y1++){
                            for (int y2 = y1; y2 <= 2; y2++) check_report(points, tree, x1, y1, x2, y2);
                        }
                    }
                }
                for (int qx = -1; qx <= 3; qx++){
                    for (int qy = -1; qy <= 2; qy++){
                        for (int skip = -1; skip < n; skip++) check_nearest(points, tree, {qx, qy}, skip);
                    }
                }
            }
        }
    }
}

/// Points at corners and edges of the |x|, |y| <= 1e9 box, where a squared distance reaches 8e18
vector<Point> extreme_points(int n){
    vector<Point> points;
    vector<long long> coords = {-C, -C + 1, 0, C - 1, C};
    for (int i = 0; i < n; i++){
        points.push_back({coords[stress::rand_int(0, 4)], coords[stress::rand_int(0, 4)]});
    }

    return points;
}

/// Points on a circle around the origin, the worst case: a query at the centre has every point near the best
vector<Point> ring_points(int n){
    vector<Point> points;
    for (int i = 0; i < n; i++){
        double angle = 2 * acos(-1.0) * i / n;
        points.push_back({llround(1e8 * cos(angle)), llround(1e8 * sin(angle))});
    }

    return points;
}

vector<Point> clustered_points(int n){
    vector<Point> centres;
    for (int k = stress::rand_int(1, 5); k > 0; k--){
        centres.push_back({stress::rand_int(-C / 2, C / 2), stress::rand_int(-C / 2, C / 2)});
    }

    vector<Point> points;
    for (int i = 0; i < n; i++){
        Point c = centres[stress::rand_int(0, centres.size() - 1)];
        long long spread = vector<long long>{0, 3, 1000}[stress::rand_int(0, 2)];
        points.push_back({c.x + stress::rand_int(-spread, spread), c.y + stress::rand_int(-spread, spread)});
    }

    return points;
}

/// Compares nearest(q, skip) and report with O(n) scans, both metrics, on random, degenerate and adversarial point sets
int main(){
    check_exhaustive();

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 50 ? 60 : 1000);
        long long range = vector<long long>{0, 2, 10, 1000, C}[stress::rand_int(0, 4)];

        vector<Point> points;
        if (it % 7 == 0) points = extreme_points(n);
        else if (it % 7 == 1) points = clustered_points(n);
        else{
            long long fixed = stress::rand_int(-range, range);
            for (int i = 0; i < n; i++){
                Point p = {stress::rand_int(-range, range), stress::rand_int(-range, range)};
                if (it % 7 == 2) p.x = fixed;
                if (it % 7 == 3) p.y = fixed;
                if (it % 7 == 4) p.y = p.x;
                points.push_back(p);
            }
        }
        check_all(points, 20, range == 0 ? 3 : (it % 7 == 0 ? C : range));
    }

    vector<Point> grid;
    for (int x = 0; x < 40; x++){
        for (int y = 0; y < 40; y++) grid.push_back({x * 7, y * 7});
    }
    check_all(grid, 200, 300);
    check_all(ring_points(2000), 200, 100000000);
    check_all(vector<Point>(1000, Point{-C, C}), 200, C);

    vector<Point> large;
    for (int i = 0; i < 200000; i++) large.push_back({stress::rand_int(-C, C), stress::rand_int(-C, C)});
    KDTree large_tree(large);
    for (int k = 0; k < 100; k++){
        int i = stress::rand_int(0, large.size() - 1);
        check_nearest(large, large_tree, large[i], i);
        long long x = stress::rand_int(-C, C - 50000000), y = stress::rand_int(-C, C - 50000000);
        check_report(large, large_tree, x, y, x + 50000000, y + stress::rand_int(0, 50000000));
    }

    return 0;
}
