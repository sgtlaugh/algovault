#include "../common.h"

#define main library_main
#include "../../code_library/geometry/segment_intersection_sweep.cpp"
#undef main

using i128 = __int128;

i128 cross2(long long ax, long long ay, long long bx, long long by){
    return (i128)ax * by - (i128)ay * bx;
}

/// Solves p.a + t * r = q.a + u * s for t, u in [0, 1] as exact fractions, collinear pairs by overlap of projections
bool brute_intersect(const Segment& p, const Segment& q){
    long long rx = p.b.x - p.a.x, ry = p.b.y - p.a.y, sx = q.b.x - q.a.x, sy = q.b.y - q.a.y;
    long long wx = q.a.x - p.a.x, wy = q.a.y - p.a.y;
    i128 den = cross2(rx, ry, sx, sy);

    if (den != 0){
        i128 t = cross2(wx, wy, sx, sy), u = cross2(wx, wy, rx, ry);
        if (den < 0) den = -den, t = -t, u = -u;
        return 0 <= t && t <= den && 0 <= u && u <= den;
    }

    if (cross2(wx, wy, rx, ry) != 0 || cross2(wx, wy, sx, sy) != 0) return false;
    if (rx == 0 && ry == 0 && sx == 0 && sy == 0) return wx == 0 && wy == 0;

    long long dx = (rx || ry) ? rx : sx, dy = (rx || ry) ? ry : sy;
    auto proj = [&](const Point& v){ return (i128)(v.x - p.a.x) * dx + (i128)(v.y - p.a.y) * dy; };
    i128 p_lo = min(proj(p.a), proj(p.b)), p_hi = max(proj(p.a), proj(p.b));
    i128 q_lo = min(proj(q.a), proj(q.b)), q_hi = max(proj(q.a), proj(q.b));
    return max(p_lo, q_lo) <= min(p_hi, q_hi);
}

bool brute_any(const vector<Segment>& segs){
    for (int i = 0; i < (int)segs.size(); i++){
        for (int j = i + 1; j < (int)segs.size(); j++){
            if (brute_intersect(segs[i], segs[j])) return true;
        }
    }
    return false;
}

void check_pair(const vector<Segment>& segs, pair<int, int> res, bool expected){
    if (!expected){
        assert(res == make_pair(-1, -1));
        return;
    }
    assert(0 <= res.first && res.first < res.second && res.second < (int)segs.size());
    assert(brute_intersect(segs[res.first], segs[res.second]));
}

void check(const vector<Segment>& segs){
    check_pair(segs, find_intersecting_pair(segs), brute_any(segs));
}

Point rand_point(long long lo, long long hi){
    return {stress::rand_int(lo, hi), stress::rand_int(lo, hi)};
}

/// Mixes points, vertical, horizontal and short and long slanted segments, all inside [lo, hi]^2
Segment rand_segment(long long lo, long long hi){
    Point a = rand_point(lo, hi);
    long long len = stress::rand_int(0, 2) ? max(1LL, (hi - lo) / 8) : hi - lo;
    auto near = [&](long long v){ return min(hi, max(lo, v + stress::rand_int(-len, len))); };

    int kind = stress::rand_int(0, 5);
    if (kind == 0) return {a, a};
    if (kind == 1) return {a, {a.x, near(a.y)}};
    if (kind == 2) return {a, {near(a.x), a.y}};
    return {a, {near(a.x), near(a.y)}};
}

/// Random segments kept only when they miss every segment kept so far
vector<Segment> disjoint_set(int target, long long lo, long long hi){
    vector<Segment> segs;
    for (int attempt = 0; attempt < 6 * target && (int)segs.size() < target; attempt++){
        Segment cand = rand_segment(lo, hi);
        bool ok = true;
        for (auto& s : segs) if (brute_intersect(s, cand)){ ok = false; break; }
        if (ok) segs.push_back(cand);
    }
    return segs;
}

vector<Segment> scaled_up(vector<Segment> segs, long long k){
    for (auto& s : segs) for (Point* p : {&s.a, &s.b}) p->x *= k, p->y *= k;
    return segs;
}

void check_variants(vector<Segment> segs, long long lo, long long hi){
    shuffle(segs.begin(), segs.end(), stress::rng());
    for (auto& s : segs) if (stress::rand_int(0, 1)) swap(s.a, s.b);
    check(segs);

    if (segs.empty()) return;
    vector<Segment> dup = segs;
    dup.insert(dup.begin() + stress::rand_int(0, dup.size()), segs[stress::rand_int(0, segs.size() - 1)]);
    check_pair(dup, find_intersecting_pair(dup), true);

    vector<Segment> extra = segs;
    extra.insert(extra.begin() + stress::rand_int(0, extra.size()), rand_segment(lo, hi));
    check(extra);
}

void exhaustive_small(){
    vector<Point> grid9, grid6;
    for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) grid9.push_back({x, y});
    for (int x = 0; x < 3; x++) for (int y = 0; y < 2; y++) grid6.push_back({x, y});

    vector<Segment> segs9, segs6;
    for (auto& a : grid9) for (auto& b : grid9) segs9.push_back({a, b});
    for (auto& a : grid6) for (auto& b : grid6) segs6.push_back({a, b});

    for (auto& s : segs9) for (auto& t : segs9) check({s, t});
    for (auto& s : segs6) for (auto& t : segs6) for (auto& u : segs6) check({s, t, u});
}

/// n segments on distinct y bands never intersect, then one vertical segment crosses a known range of bands
void large_bands(int n, long long y0, long long gap, long long x_range){
    vector<Segment> segs;
    for (int i = 0; i < n; i++){
        long long x1 = stress::rand_int(-x_range, x_range), x2 = stress::rand_int(-x_range, x_range), y = y0 + i * gap;
        segs.push_back({{x1, y}, {x2, y + stress::rand_int(0, gap / 2)}});
    }
    shuffle(segs.begin(), segs.end(), stress::rng());
    assert(find_intersecting_pair(segs) == make_pair(-1, -1));

    Segment pole = {{segs[0].a.x, y0}, {segs[0].a.x, y0 + (n - 1) * gap}};
    segs.insert(segs.begin() + stress::rand_int(0, n), pole);
    check_pair(segs, find_intersecting_pair(segs), true);
}

int main(){
    const long long B = 1000000000000000000LL, K = B / 4;

    exhaustive_small();

    for (long long it = 0; it < stress::scaled(30000); it++){
        int n = stress::rand_int(0, 8);
        long long range = it % 3 == 0 ? 2 : it % 3 == 1 ? 4 : 10;
        vector<Segment> segs;
        for (int i = 0; i < n; i++) segs.push_back(rand_segment(0, range));
        check(segs);
        if (it % 4 == 0) check(scaled_up(segs, B / range));
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int target = it % 500 == 0 ? 2000 : it % 10 == 0 ? 400 : stress::rand_int(1, 60);
        long long hi = it % 3 == 0 ? 6 : it % 3 == 1 ? 30 : 1000000;
        check_variants(disjoint_set(target, 0, hi), 0, hi);
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        vector<Segment> segs = disjoint_set(stress::rand_int(1, 40), -4, 4);
        check_variants(scaled_up(segs, K), -B, B);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        vector<Segment> segs;
        for (int i = 0, n = stress::rand_int(2, 6); i < n; i++) segs.push_back(rand_segment(-B, B));
        check(segs);
    }

    large_bands(200000, -1000000000, 3, 1000000000);
    large_bands(200000, -B, B / 100000, B);

    return 0;
}
