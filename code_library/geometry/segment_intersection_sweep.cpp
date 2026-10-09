/***
 *
 * Segment Intersection Sweep (Shamos-Hoey)
 * Decides whether any two of n closed segments intersect and returns one intersecting pair
 *
 * Complexity: O(n log n), O(n) memory
 *
 * find_intersecting_pair(segs) returns {i, j} with i < j such that segs[i] and segs[j] share a point, or {-1, -1}
 *   Endpoints are inclusive: touching, collinear overlapping, duplicate and single-point segments all intersect
 *   Exact integer predicates in __int128, |x|, |y| <= 1e18
 *   Which pair is returned when several intersect is unspecified
 *
 * Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;

    bool operator<(const Point& p) const{
        return x < p.x || (x == p.x && y < p.y);
    }
};

struct Segment{
    Point a, b;
};

int orientation(const Point& o, const Point& a, const Point& b){
    __int128 v = (__int128)(a.x - o.x) * (b.y - o.y) - (__int128)(a.y - o.y) * (b.x - o.x);
    return (v > 0) - (v < 0);
}

bool on_segment(const Point& p, const Point& a, const Point& b){
    return orientation(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

bool segments_intersect(const Segment& s, const Segment& t){
    int d1 = orientation(s.a, s.b, t.a), d2 = orientation(s.a, s.b, t.b), d3 = orientation(t.a, t.b, s.a), d4 = orientation(t.a, t.b, s.b);
    if (d1 * d2 < 0 && d3 * d4 < 0) return true;
    return on_segment(t.a, s.a, s.b) || on_segment(t.b, s.a, s.b) || on_segment(s.a, t.a, t.b) || on_segment(s.b, t.a, t.b);
}

pair<int, int> find_intersecting_pair(const vector<Segment>& segs){
    int n = segs.size();
    vector<Segment> s = segs;
    vector<tuple<long long, long long, int, int>> events;
    events.reserve(2 * n);

    /// Sweeping points in (x, y) order tilts the sweep line infinitesimally, so vertical segments need no special case
    /// Insertions (type 0) precede removals at the same point, or segments touching only there would be missed
    for (int i = 0; i < n; i++){
        if (s[i].b < s[i].a) swap(s[i].a, s[i].b);
        events.emplace_back(s[i].a.x, s[i].a.y, 0, i);
        events.emplace_back(s[i].b.x, s[i].b.y, 1, i);
    }
    sort(events.begin(), events.end());

    /// Bottom to top along the sweep line through the later left endpoint, it is active in both
    /// Neither is below the other only when they share a point, and lower_bound then hands that segment to the check
    auto below = [&](int i, int j){
        if (s[i].a < s[j].a) return orientation(s[i].a, s[i].b, s[j].a) > 0;
        return orientation(s[j].a, s[j].b, s[i].a) < 0;
    };
    set<int, decltype(below)> active(below);
    vector<set<int, decltype(below)>::iterator> where(n);
    auto ordered = [](int i, int j){ return make_pair(min(i, j), max(i, j)); };

    for (auto& e : events){
        int type = get<2>(e), id = get<3>(e);
        if (type == 0){
            auto above = active.lower_bound(id);
            if (above != active.end() && segments_intersect(s[*above], s[id])) return ordered(*above, id);
            if (above != active.begin() && segments_intersect(s[*prev(above)], s[id])) return ordered(*prev(above), id);
            where[id] = active.insert(above, id);
            continue;
        }

        auto it = where[id];
        if (it != active.begin() && next(it) != active.end() && segments_intersect(s[*prev(it)], s[*next(it)])){
            return ordered(*prev(it), *next(it));
        }
        active.erase(it);
    }

    return {-1, -1};
}

int main(){
    const long long B = 1000000000000000000LL;
    auto none = make_pair(-1, -1), first_two = make_pair(0, 1);

    assert(find_intersecting_pair({}) == none);
    assert(find_intersecting_pair({{{0, 0}, {5, 5}}}) == none);
    assert(find_intersecting_pair({{{0, 0}, {4, 4}}, {{0, 4}, {4, 0}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {2, 2}}, {{5, 0}, {2, 2}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {4, 0}}, {{2, 3}, {2, 0}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {4, 0}}, {{0, 1}, {4, 1}}}) == none);
    assert(find_intersecting_pair({{{0, 0}, {1, 1}}, {{3, 3}, {2, 2}}}) == none);
    assert(find_intersecting_pair({{{0, 0}, {2, 2}}, {{1, 1}, {3, 3}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {0, 1}}, {{0, 3}, {0, 2}}}) == none);
    assert(find_intersecting_pair({{{0, 0}, {0, 2}}, {{0, 2}, {0, 3}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {2, 2}}, {{1, 1}, {1, 1}}}) == first_two);
    assert(find_intersecting_pair({{{0, 0}, {2, 2}}, {{1, 2}, {1, 2}}}) == none);
    assert(find_intersecting_pair({{{7, 7}, {7, 7}}, {{7, 7}, {7, 7}}}) == first_two);

    assert(find_intersecting_pair({{{0, 0}, {10, 0}}, {{0, 2}, {10, 2}}, {{0, 4}, {10, 4}}, {{3, 1}, {5, 3}}}) == make_pair(1, 3));
    assert(find_intersecting_pair({{{0, 0}, {10, 10}}, {{1, 10}, {10, 0}}, {{0, 5}, {2, 5}}}) == first_two);
    assert(find_intersecting_pair({{{0, 5}, {2, 5}}, {{0, 0}, {10, 10}}, {{6, 9}, {7, 9}}, {{1, 10}, {10, 0}}}) == make_pair(1, 3));

    assert(find_intersecting_pair({{{-B, -B}, {B, B}}, {{-B, B}, {B, -B}}}) == first_two);
    assert(find_intersecting_pair({{{-B, -B}, {B - 1, B - 1}}, {{-B, -B + 1}, {B - 1, B}}}) == none);
    assert(find_intersecting_pair({{{-B, 0}, {B, 0}}, {{0, 1}, {B, B}}}) == none);
    assert(find_intersecting_pair({{{-B, 0}, {B, 0}}, {{0, 0}, {B, B}}}) == first_two);
    assert(find_intersecting_pair({{{-B, -B}, {B, B - 1}}, {{B - 2, B - 2}, {B - 2, B - 2}}}) == none);

    return 0;
}
