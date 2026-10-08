#include "common.h"

#define main library_main
#include "../code_library/convex_hull.cpp"
#undef main

bool same(const Point& a, const Point& b){
    return a.x == b.x && a.y == b.y;
}

/// A strict hull is the unique strictly convex counter-clockwise polygon over input points that contains every input point
void check_hull(const vector<Point>& points, const vector<Point>& hull){
    vector<Point> distinct = points;
    sort(distinct.begin(), distinct.end());
    distinct.erase(unique(distinct.begin(), distinct.end(), same), distinct.end());
    int h = hull.size();

    if (distinct.size() <= 1){
        assert(h == (int)distinct.size() && (h == 0 || same(hull[0], distinct[0])));
        return;
    }
    assert(h >= 2 && same(hull[0], distinct[0]));
    for (auto& v : hull) assert(binary_search(distinct.begin(), distinct.end(), v));
    for (int i = 0; i < h; i++) for (int j = i + 1; j < h; j++) assert(!same(hull[i], hull[j]));

    if (h == 2){  /// every point lies on the segment between the two endpoints
        assert(same(hull[1], distinct.back()));
        for (auto& p : points) assert(cross(hull[0], hull[1], p) == 0);
        return;
    }
    for (int i = 0; i < h; i++){
        assert(cross(hull[i], hull[(i + 1) % h], hull[(i + 2) % h]) > 0);
        for (auto& p : points) assert(cross(hull[i], hull[(i + 1) % h], p) >= 0);
    }
}

int main(){
    /// Folds along an axis-parallel edge, where only one coordinate changes direction
    assert(!is_convex({Point(0, 0), Point(3, 0), Point(1, 0), Point(4, 0), Point(0, 4)}));
    assert(!is_convex({Point(0, 0), Point(4, 0), Point(0, 4), Point(0, 1), Point(0, 3)}));

    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 300), range = stress::rand_int(0, 2) ? 4 : 1000000000;
        vector<Point> points;
        for (int i = 0; i < n; i++){
            if (stress::rand_int(0, 4) == 0 && i) points.push_back(points[stress::rand_int(0, i - 1)]);              /// duplicates
            else if (stress::rand_int(0, 5) == 0) points.push_back(Point(stress::rand_int(-range, range), 7));    /// collinear row
            else points.push_back(Point(stress::rand_int(-range, range), stress::rand_int(-range, range)));
        }

        auto hull = get_convex_hull(points);
        check_hull(points, hull);

        int h = hull.size();
        if (h < 3){
            assert(!is_convex(hull));
            continue;
        }

        vector<Point> rotated(h);
        int shift = stress::rand_int(0, h - 1);
        for (int i = 0; i < h; i++) rotated[i] = hull[(i + shift) % h];
        assert(is_convex(rotated) && is_convex(vector<Point>(rotated.rbegin(), rotated.rend())));

        /// Inserting the midpoint of an edge keeps it convex, when the midpoint is a lattice point
        Point a = hull[0], b = hull[1];
        if ((a.x + b.x) % 2 == 0 && (a.y + b.y) % 2 == 0){
            vector<Point> with_mid = hull;
            with_mid.insert(with_mid.begin() + 1, Point((a.x + b.x) / 2, (a.y + b.y) / 2));
            assert(is_convex(with_mid));
        }

        if (range == 4){  /// an edge doubling back on itself, a -> 2/3 -> 1/3 -> b, turns one way but is not a simple polygon
            vector<Point> folded;
            for (auto& p : hull) folded.push_back(Point(3 * p.x, 3 * p.y));
            Point s = folded[0], t = folded[1], d((t.x - s.x) / 3, (t.y - s.y) / 3);
            folded.insert(folded.begin() + 1, {Point(s.x + 2 * d.x, s.y + 2 * d.y), Point(s.x + d.x, s.y + d.y)});
            assert(!is_convex(folded));
        }

        if (h >= 4){  /// swapping two adjacent vertices makes the boundary cross itself
            vector<Point> crossed = hull;
            int i = stress::rand_int(0, h - 1);
            swap(crossed[i], crossed[(i + 1) % h]);
            assert(!is_convex(crossed));
        }
        if (h >= 5){  /// visiting every k-th vertex for k coprime with h winds around more than once: a star
            for (int k = 2; k < h - 1; k++){
                if (__gcd(k, h) != 1) continue;
                vector<Point> star;
                for (int i = 0; i < h; i++) star.push_back(hull[(long long)i * k % h]);
                assert(!is_convex(star));
            }
        }
    }
    return 0;
}
