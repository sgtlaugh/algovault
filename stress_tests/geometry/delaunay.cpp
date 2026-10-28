#include "../common.h"

#define main library_main
#include "../../code_library/geometry/delaunay.cpp"
#undef main

const long long LIMIT = 1000000000;

__int128 cross128(const Point& o, const Point& a, const Point& b){
    return (__int128)(a.x - o.x) * (b.y - o.y) - (__int128)(a.y - o.y) * (b.x - o.x);
}

long long dist2(const Point& a, const Point& b){
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

/// Sign of d against the circle through ccw a, b, c, the lifted determinant taken around a instead of d
int circle_side(const Point& a, const Point& b, const Point& c, const Point& d){
    __int128 lb = dist2(b, a), lc = dist2(c, a), ld = dist2(d, a);
    __int128 det = lb * cross128(a, c, d) + lc * cross128(a, d, b) + ld * cross128(a, b, c);
    return (det < 0) - (det > 0);
}

bool before(const Point& a, const Point& b){
    return make_pair(a.x, a.y) < make_pair(b.x, b.y);
}

/// Points of the strict hull, ccw, via gift wrapping in 128-bit arithmetic
vector<Point> brute_hull(vector<Point> p){
    sort(p.begin(), p.end(), before);
    vector<Point> res;
    Point cur = p[0];
    do {
        res.push_back(cur);
        Point next = (p[0].x == cur.x && p[0].y == cur.y) ? p[1] : p[0];
        for (auto& q : p){
            __int128 c = cross128(cur, next, q);
            if (c < 0 || (c == 0 && dist2(cur, q) > dist2(cur, next))) next = q;
        }
        cur = next;
    } while (cur.x != p[0].x || cur.y != p[0].y);
    return res;
}

bool on_segment(const Point& p, const Point& a, const Point& b){
    return cross128(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x)
        && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

/// Checks everything the header promises; full_circle_check tests every point against every circumcircle
void validate(const vector<Point>& points, const Delaunay& dt, bool full_circle_check){
    int n = points.size();
    map<pair<long long, long long>, int> first;
    for (int i = 0; i < n; i++) first.insert({{points[i].x, points[i].y}, i});
    vector<int> rep(n), uniq;
    for (int i = 0; i < n; i++){
        rep[i] = first[{points[i].x, points[i].y}];
        if (rep[i] == i) uniq.push_back(i);
    }
    int m = uniq.size();

    set<pair<int, int>> base;
    int duplicates = 0;
    for (auto [u, v] : dt.edges){
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        if (rep[u] != u){
            assert(v == rep[u]);  /// a duplicate joined to anything but its smallest-index copy
            duplicates++;
            continue;
        }
        assert(rep[v] == v);
        assert(base.insert(minmax(u, v)).second);
    }
    assert(duplicates == n - m);

    bool collinear = true;
    for (int i = 2; i < m; i++) collinear &= cross128(points[uniq[0]], points[uniq[1]], points[uniq[i]]) == 0;
    if (collinear){
        assert(dt.triangles.empty());
        vector<int> line = uniq;
        sort(line.begin(), line.end(), [&](int a, int b){ return before(points[a], points[b]); });
        set<pair<int, int>> chain;
        for (int i = 0; i + 1 < m; i++) chain.insert(minmax(line[i], line[i + 1]));
        assert(base == chain);
        return;
    }

    set<pair<int, int>> directed;
    vector<bool> used(n, false);
    __int128 area = 0;
    for (auto [a, b, c] : dt.triangles){
        assert(rep[a] == a && rep[b] == b && rep[c] == c);
        assert(a < b && a < c);
        __int128 turn = cross128(points[a], points[b], points[c]);
        assert(turn > 0);
        area += turn;
        used[a] = used[b] = used[c] = true;
        for (auto [u, v] : {make_pair(a, b), make_pair(b, c), make_pair(c, a)}){
            assert(directed.insert({u, v}).second);  /// two triangles overlapping across an edge
        }
    }
    for (int i : uniq) assert(used[i]);

    set<pair<int, int>> from_triangles;
    int boundary = 0;
    for (auto [u, v] : directed){
        from_triangles.insert(minmax(u, v));
        if (directed.count({v, u})) continue;
        boundary++;
        for (int i : uniq) assert(cross128(points[u], points[v], points[i]) >= 0);
    }
    assert(from_triangles == base);

    vector<Point> distinct;
    for (int i : uniq) distinct.push_back(points[i]);
    vector<Point> hull = brute_hull(distinct);
    int h = hull.size();
    __int128 hull_area = 0;
    for (int i = 0; i < h; i++) hull_area += cross128(hull[0], hull[i], hull[(i + 1) % h]);
    assert(area == hull_area);

    int on_hull = 0;
    for (auto& p : distinct){
        bool on = false;
        for (int i = 0; i < h && !on; i++) on = on_segment(p, hull[i], hull[(i + 1) % h]);
        on_hull += on;
    }
    assert(boundary == on_hull);
    assert((int)dt.triangles.size() == 2 * m - 2 - on_hull);

    if (full_circle_check){
        for (auto [a, b, c] : dt.triangles){
            for (int i : uniq) assert(circle_side(points[a], points[b], points[c], points[i]) <= 0);
        }

        return;
    }

    /// Locally Delaunay at every interior edge implies globally Delaunay
    map<pair<int, int>, int> apex;
    for (auto [a, b, c] : dt.triangles){
        apex[{a, b}] = c, apex[{b, c}] = a, apex[{c, a}] = b;
    }
    for (auto [e, c] : apex){
        auto it = apex.find({e.second, e.first});
        if (it == apex.end()) continue;
        assert(circle_side(points[e.first], points[e.second], points[c], points[it->second]) <= 0);
    }
}

/// In general position the Delaunay triangles are exactly the triples with an empty circumcircle
set<array<int, 3>> brute_triangles(const vector<Point>& p){
    int n = p.size();
    set<array<int, 3>> res;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int k = j + 1; k < n; k++){
        __int128 turn = cross128(p[i], p[j], p[k]);
        array<int, 3> t = turn > 0 ? array<int, 3>{i, j, k} : array<int, 3>{i, k, j};
        bool empty = true;
        for (int q = 0; q < n && empty; q++) empty = circle_side(p[t[0]], p[t[1]], p[t[2]], p[q]) <= 0;
        if (empty) res.insert(t);
    }

    return res;
}

bool general_position(const vector<Point>& p){
    int n = p.size();
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int k = j + 1; k < n; k++){
        __int128 turn = cross128(p[i], p[j], p[k]);
        if (turn == 0) return false;
        for (int q = k + 1; q < n; q++){
            if (turn > 0 ? circle_side(p[i], p[j], p[k], p[q]) == 0 : circle_side(p[i], p[k], p[j], p[q]) == 0){
                return false;
            }
        }
    }

    return true;
}

vector<long long> kruskal_weights(const vector<Point>& p, const vector<array<int, 2>>& edges){
    vector<array<int, 2>> sorted = edges;
    sort(sorted.begin(), sorted.end(), [&](const array<int, 2>& a, const array<int, 2>& b){
        return dist2(p[a[0]], p[a[1]]) < dist2(p[b[0]], p[b[1]]);
    });
    vector<int> parent(p.size());
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    vector<long long> res;
    for (auto [u, v] : sorted){
        if (find(u) == find(v)) continue;
        parent[find(u)] = find(v);
        res.push_back(dist2(p[u], p[v]));
    }

    return res;
}

/// The sorted multiset of MST edge weights is the same for every MST
vector<long long> prim_weights(const vector<Point>& p){
    int n = p.size();
    vector<long long> best(n, LLONG_MAX), res;
    vector<bool> done(n, false);
    best[0] = 0;
    for (int it = 0; it < n; it++){
        int u = -1;
        for (int i = 0; i < n; i++) if (!done[i] && (u < 0 || best[i] < best[u])) u = i;
        done[u] = true;
        if (it) res.push_back(best[u]);
        for (int i = 0; i < n; i++) if (!done[i]) best[i] = min(best[i], dist2(p[u], p[i]));
    }
    sort(res.begin(), res.end());
    return res;
}

void check(const vector<Point>& points, bool full_circle_check = true){
    Delaunay dt(points);
    validate(points, dt, full_circle_check);
    if (full_circle_check && !points.empty()) assert(kruskal_weights(points, dt.edges) == prim_weights(points));
}

vector<Point> random_points(int n, long long lo, long long hi){
    vector<Point> p(n);
    for (auto& q : p) q = {stress::rand_int(lo, hi), stress::rand_int(lo, hi)};
    return p;
}

/// A k by k grid with step s centered so its corners reach the coordinate limit, every unit square cocircular
vector<Point> grid(int k, long long step, double keep){
    vector<Point> p;
    long long start = -step * (k - 1) / 2;
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++){
        if (stress::rand_int(0, 999) < keep * 1000) p.push_back({start + i * step, start + j * step});
    }
    shuffle(p.begin(), p.end(), stress::rng());
    return p;
}

void add_duplicates(vector<Point>& p, int count){
    for (int i = 0; i < count && !p.empty(); i++) p.push_back(p[stress::rand_int(0, p.size() - 1)]);
    shuffle(p.begin(), p.end(), stress::rng());
}

int main(){
    check({});
    check({{5, 5}});
    check({{5, 5}, {5, 5}, {5, 5}});
    check({{-LIMIT, -LIMIT}, {LIMIT, LIMIT}});
    check({{-LIMIT, -LIMIT}, {LIMIT, -LIMIT}, {LIMIT, LIMIT}, {-LIMIT, LIMIT}, {0, 0}, {1, 0}});

    /// Every subset of a 3 by 3 grid, the smallest home of collinear and cocircular cases
    for (int mask = 0; mask < (1 << 9); mask++){
        vector<Point> p;
        for (int i = 0; i < 9; i++) if (mask >> i & 1) p.push_back({i / 3, i % 3});
        check(p);
        if (!p.empty()) add_duplicates(p, 2);
        check(p);
    }

    for (int it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(3, 8);
        long long range = stress::rand_int(0, 1) ? 30 : LIMIT;
        vector<Point> p = random_points(n, -range, range);
        if (!general_position(p)) continue;
        Delaunay dt(p);
        set<array<int, 3>> got(dt.triangles.begin(), dt.triangles.end());
        assert(got.size() == dt.triangles.size());
        assert(got == brute_triangles(p));
    }

    for (int it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 60), shape = stress::rand_int(0, 5);
        vector<Point> p;
        if (shape == 0) p = random_points(n, -3, 3);
        else if (shape == 1) p = random_points(n, -LIMIT, LIMIT);
        else if (shape == 2) p = random_points(n, LIMIT - 20, LIMIT);
        else if (shape == 3) p = grid(stress::rand_int(2, 9), 2 * LIMIT / 8, 0.7);
        else if (shape == 4){
            long long dx = stress::rand_int(-3, 3), dy = stress::rand_int(-3, 3);
            for (int i = 0; i < n; i++){
                long long t = stress::rand_int(-20, 20);
                p.push_back({t * dx, t * dy});
            }
        }
        else {
            /// Two far clusters at opposite corners, merged across a huge gap
            for (int i = 0; i < n; i++){
                long long s = i % 2 ? 1 : -1;
                p.push_back({s * stress::rand_int(LIMIT - 5, LIMIT), s * stress::rand_int(LIMIT - 5, LIMIT)});
            }
        }
        if (stress::rand_int(0, 2) == 0) add_duplicates(p, stress::rand_int(1, 5));
        check(p);
    }

    vector<vector<Point>> large = {
        random_points(stress::scaled(15000), -LIMIT, LIMIT),
        grid(81, 2 * LIMIT / 80, 1.0),
        grid(101, 2 * LIMIT / 100, 0.5),
    };
    vector<Point> line;
    for (int i = 0; i < stress::scaled(10000); i++){
        long long t = stress::rand_int(-LIMIT / 3, LIMIT / 3);
        line.push_back({t, 3 * t - 1});
    }
    large.push_back(line);
    add_duplicates(large[0], 1000);
    for (auto& p : large) check(p, false);

    return 0;
}
