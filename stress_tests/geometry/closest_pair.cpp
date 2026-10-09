#include "../common.h"

#define main library_main
#include "../../code_library/geometry/closest_pair.cpp"
#undef main

const long long C = 1000000000;

long long dist2(const Point& a, const Point& b){
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

long long brute(const vector<Point>& points){
    long long best = LLONG_MAX;
    for (int i = 0; i < (int)points.size(); i++){
        for (int j = i + 1; j < (int)points.size(); j++) best = min(best, dist2(points[i], points[j]));
    }
    return best;
}

void check(const vector<Point>& points, long long expected){
    ClosestPairResult res = closest_pair(points);
    assert(res.dist2 == expected);
    assert(0 <= res.i && res.i < res.j && res.j < (int)points.size());
    assert(dist2(points[res.i], points[res.j]) == res.dist2);
}

/// Distinct coordinates along one axis, the answer is the smallest gap between sorted neighbours
void check_collinear(int n, bool vertical){
    vector<long long> coords;
    for (int i = 0; i < n; i++) coords.push_back(stress::rand_int(-C, C));
    sort(coords.begin(), coords.end());
    coords.erase(unique(coords.begin(), coords.end()), coords.end());

    long long gap = LLONG_MAX;
    for (int i = 1; i < (int)coords.size(); i++) gap = min(gap, coords[i] - coords[i - 1]);
    shuffle(coords.begin(), coords.end(), stress::rng());

    long long other = stress::rand_int(-C, C);
    vector<Point> points;
    for (long long c : coords) points.push_back(vertical ? Point{other, c} : Point{c, other});
    check(points, gap * gap);
}

/// A full lattice with spacing 1000 (keeps the sweep box full) plus one point planted at offset (dx, dy) from a lattice point
void check_lattice(int side){
    long long dx = stress::rand_int(-20, 20), dy = stress::rand_int(1, 20);
    vector<Point> points;
    for (int i = 0; i < side; i++){
        for (int j = 0; j < side; j++) points.push_back({i * 1000LL - C, j * 1000LL - C});
    }
    Point anchor = points[stress::rand_int(0, side * side - 1)];
    points.push_back({anchor.x + dx, anchor.y + dy});
    shuffle(points.begin(), points.end(), stress::rng());

    check(points, dx * dx + dy * dy);
}

/// Compares the squared distance against an O(n^2) scan and checks the returned indices form a pair at that distance
int main(){
    for (long long it = 0; it < stress::scaled(30000); it++){
        int n = it < 3000 ? 2 + it % 7 : stress::rand_int(2, it % 50 ? 60 : 1500);
        long long range = vector<long long>{2, 10, 1000, C}[stress::rand_int(0, 3)];

        vector<Point> points;
        for (int i = 0; i < n; i++){
            if (it % 4 == 0 && stress::rand_int(0, 9) == 0) points.push_back({stress::rand_int(0, 1) ? -C : C, stress::rand_int(0, 1) ? -C : C});
            else if (stress::rand_int(0, 9) == 0) points.push_back({stress::rand_int(-range, range), 3});
            else points.push_back({stress::rand_int(-range, range), stress::rand_int(-range, range)});
        }
        check(points, brute(points));
    }

    for (int it = 0; it < 4; it++){
        check_collinear(200000, it % 2);
        check_lattice(300);
    }

    return 0;
}
