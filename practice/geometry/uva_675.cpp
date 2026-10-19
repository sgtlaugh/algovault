/***
 *
 * Convex Hull
 * Monotone chain hull, rotating calipers, Minkowski sum, O(log n) queries on a hull and an online insertion hull
 *
 * Complexity: O(n log n) for get_convex_hull, the rest as listed below
 *
 * Coordinates satisfy |x|, |y| <= 1e9, so every cross product, dot product and squared distance fits in int64_t
 * A "hull" below is a strictly convex polygon in counter-clockwise order, as returned by get_convex_hull
 *
 *   get_convex_hull(points): the strict hull, O(n log n)
 *   is_convex(polygon): O(n)
 *   farthest_pair(hull): two vertices at maximum distance, dist2 of them is the exact squared diameter, O(n), n >= 1
 *   width(hull): minimum distance between two parallel lines enclosing the hull, O(n), 0 when n <= 2
 *   minkowski_sum(a, b): the strict hull of {p + q}, O(n + m), a and b convex counter-clockwise polygons
 *       with any starting vertex, collinear vertices allowed, result coordinates reach 2e9
 *   extreme_vertex(hull, dir): index of a vertex maximizing dot(vertex, dir), O(log n), n >= 1, |dir.x|, |dir.y| <= 2e9
 *   tangents(hull, q): {i, j}, the hull lies left of the ray q -> hull[i] and right of q -> hull[j], O(log n)
 *       q strictly outside the hull, either endpoint is returned when q is collinear with an edge
 *   line_hull(a, b, hull): the line through a != b against a hull with n >= 3, O(log n)
 *       {-1, -1} no intersection, {i, -1} touches only vertex i, {i, i} contains edge (i, i + 1)
 *       {i, j} crosses edges (i, i + 1) then (j, j + 1) in the direction a -> b, a crossed vertex k counts as edge (k, k + 1)
 *   DynamicHull: add(p) amortized O(log n), contains(p) O(log n) boundary inclusive, hull() O(h) as get_convex_hull returns it
 *
 * Optimization notes: get_convex_hull can be converted to O(n) using radix sort
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    int64_t x, y;

    Point() {}
    Point(int64_t x, int64_t y) : x(x), y(y) {}

    inline bool operator < (const Point &p) const {
        return ((x < p.x) || (x == p.x && y < p.y));
    }

    bool operator == (const Point& p) const {
        return x == p.x && y == p.y;
    }

    Point operator + (const Point& p) const {
        return Point(x + p.x, y + p.y);
    }

    Point operator - (const Point& p) const {
        return Point(x - p.x, y - p.y);
    }
};

int sgn(int64_t v){
    return (v > 0) - (v < 0);
}

int64_t cross(const Point& a, const Point& b){
    return a.x * b.y - a.y * b.x;
}

int64_t cross(const Point &O, const Point &A, const Point &B){
    return ((A.x - O.x) * (B.y - O.y)) - ((A.y - O.y) * (B.x - O.x));
}

int64_t dot(const Point& a, const Point& b){
    return a.x * b.x + a.y * b.y;
}

int64_t dist2(const Point& a, const Point& b){
    return dot(a - b, a - b);
}

/***
 *
 * Returns the strict convex hull in counter-clockwise order, starting from the lowest x (then lowest y) point
 * Collinear and duplicate points are dropped, all-collinear input gives the two endpoints
 *
***/

vector<Point> get_convex_hull(vector<Point> P){
    sort(P.begin(), P.end());
    P.erase(unique(P.begin(), P.end()), P.end());

    int i, t, k = 0, n = P.size();
    if (n <= 1) return P;
    vector<Point> H(n << 1);

    for (i = 0; i < n; i++){
        while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }
    for (i = n - 2, t = k + 1; i >= 0; i--){
        while (k >= t && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    H.resize(k - 1);
    return H;
}

/***
 *
 * Points in P are given in clockwise or anti-clockwise order, collinear vertices are allowed
 *
***/

bool is_convex(const vector<Point>& P){
    int n = P.size(), sign = 0, first[2] = {0, 0}, last[2] = {0, 0}, flips[2] = {0, 0};
    if (n <= 2) return false; /// Line or point is not convex

    for (int i = 0; i < n; i++){
        const Point &a = P[i], &b = P[(i + 1) % n], &c = P[(i + 2) % n];
        int64_t turn = cross(a, b, c);
        if (!turn && dot(b - a, c - b) < 0) return false;
        if (turn){
            if (sign && (turn > 0) != (sign > 0)) return false;
            sign = turn > 0 ? 1 : -1;
        }

        int dir[2] = {(b.x > a.x) - (b.x < a.x), (b.y > a.y) - (b.y < a.y)};
        for (int k = 0; k < 2; k++){
            if (!dir[k]) continue;
            if (last[k] && dir[k] != last[k]) flips[k]++;
            if (!first[k]) first[k] = dir[k];
            last[k] = dir[k];
        }
    }
    for (int k = 0; k < 2; k++) flips[k] += first[k] != last[k];

    /// Same-sign turns without back-tracking still allow stars, a simple convex polygon reverses x and y exactly twice each
    return sign != 0 && flips[0] == 2 && flips[1] == 2;
}

/// Rotating calipers: calls visit(i, j) for every edge i -> i + 1 with j its first farthest vertex
template<typename Visit>
void for_each_antipode(const vector<Point>& hull, Visit visit){
    int n = hull.size();

    for (int i = 0, j = 1 % n; i < n; i++){
        const Point &a = hull[i], &b = hull[(i + 1) % n];
        while (cross(a, b, hull[(j + 1) % n]) > cross(a, b, hull[j])) j = (j + 1) % n;
        visit(i, j);
    }
}

array<Point, 2> farthest_pair(const vector<Point>& hull){
    int n = hull.size();
    array<Point, 2> best = {hull[0], hull[0]};

    for_each_antipode(hull, [&](int i, int j){
        for (int k : {i, (i + 1) % n}){
            if (dist2(hull[k], hull[j]) > dist2(best[0], best[1])) best = {hull[k], hull[j]};
        }
    });
    return best;
}

long double width(const vector<Point>& hull){
    int n = hull.size();
    if (n <= 2) return 0;
    long double best = numeric_limits<long double>::infinity();

    for_each_antipode(hull, [&](int i, int j){
        const Point &a = hull[i], &b = hull[(i + 1) % n];
        best = min(best, cross(a, b, hull[j]) / sqrtl(dist2(a, b)));
    });
    return best;
}

vector<Point> minkowski_sum(const vector<Point>& a, const vector<Point>& b){
    if (a.empty() || b.empty()) return {};

    /// Edges from the lowest vertex sort by angle in (-90, 270] degrees, so two such lists merge in O(n + m)
    auto edges = [](const vector<Point>& p){
        int n = p.size(), s = min_element(p.begin(), p.end()) - p.begin();
        vector<Point> res;
        for (int i = 0; i < n; i++){
            Point d = p[(s + i + 1) % n] - p[(s + i) % n];
            if (d.x || d.y) res.push_back(d);
        }
        return res;
    };
    auto upper_half = [](const Point& d){ return d.x > 0 || (d.x == 0 && d.y > 0); };
    auto by_angle = [&](const Point& u, const Point& v){
        if (upper_half(u) != upper_half(v)) return upper_half(u);
        return cross(u, v) > 0;
    };

    vector<Point> ea = edges(a), eb = edges(b), merged(ea.size() + eb.size());
    merge(ea.begin(), ea.end(), eb.begin(), eb.end(), merged.begin(), by_angle);

    /// Compares against the last single edge, an accumulated run of collinear edges could overflow the cross product
    vector<Point> res = {*min_element(a.begin(), a.end()) + *min_element(b.begin(), b.end())};
    Point last(0, 0);
    for (const Point& d : merged){
        if (res.size() >= 2 && cross(last, d) == 0 && dot(last, d) > 0) res.back() = res.back() + d;
        else res.push_back(res.back() + d);
        last = d;
    }

    if (res.size() > 1) res.pop_back();
    return res;
}

/// Maximum of a cyclically unimodal order on the n vertices, better(i, j) is the sign of vertex i beating vertex j
template<typename Better>
int extreme_search(int n, Better better){
    auto is_extreme = [&](int i){ return better((i + 1) % n, i) <= 0 && better(i, (i + n - 1) % n) > 0; };
    if (is_extreme(0)) return 0;
    int lo = 0, hi = n;

    while (lo + 1 < hi){
        int m = (lo + hi) / 2;
        if (is_extreme(m)) return m;
        int ls = better((lo + 1) % n, lo), ms = better((m + 1) % n, m);
        if (ls > ms || (ls == ms && ls == better(lo, m))) hi = m;
        else lo = m;
    }
    return lo;
}

int extreme_vertex(const vector<Point>& hull, const Point& dir){
    return extreme_search(hull.size(), [&](int i, int j){ return sgn(dot(hull[i] - hull[j], dir)); });
}

/// Seen from an outside point the hull spans less than 180 degrees, so the angular order behaves like a projection
pair<int, int> tangents(const vector<Point>& hull, const Point& q){
    int n = hull.size();
    int right = extreme_search(n, [&](int i, int j){ return sgn(cross(q, hull[i], hull[j])); });
    int left = extreme_search(n, [&](int i, int j){ return sgn(cross(q, hull[j], hull[i])); });
    return {right, left};
}

array<int, 2> line_hull(const Point& a, const Point& b, const vector<Point>& hull){
    int n = hull.size();
    auto side = [&](int i){ return sgn(cross(a, hull[i], b)); };
    int end_a = extreme_vertex(hull, Point(b.y - a.y, a.x - b.x));
    int end_b = extreme_vertex(hull, Point(a.y - b.y, b.x - a.x));
    if (side(end_a) < 0 || side(end_b) > 0) return {-1, -1};

    array<int, 2> res;
    for (int k = 0; k < 2; k++){
        int lo = end_b, hi = end_a;
        while ((lo + 1) % n != hi){
            int m = ((lo + hi + (lo < hi ? 0 : n)) / 2) % n;
            (side(m) == side(end_b) ? lo : hi) = m;
        }
        res[k] = (lo + !side(hi)) % n;
        swap(end_a, end_b);
    }

    if (res[0] == res[1]) return {res[0], -1};
    if (!side(res[0]) && !side(res[1])){
        int gap = (res[0] - res[1] + n + 1) % n;
        if (gap == 0) return {res[0], res[0]};
        if (gap == 2) return {res[1], res[1]};
    }
    return res;
}

struct DynamicHull{
    void add(const Point& p){
        upper.add(p);
        lower.add(Point(p.x, -p.y));
    }

    bool contains(const Point& p) const {
        return upper.covers(p) && lower.covers(Point(p.x, -p.y));
    }

    vector<Point> hull() const {
        vector<Point> res;
        for (auto& e : lower.pts) res.push_back(Point(e.first, -e.second));

        for (auto it = upper.pts.rbegin(); it != upper.pts.rend(); it++){
            Point p(it->first, it->second);
            if (!(p == res.back()) && !(p == res[0])) res.push_back(p);
        }
        return res;
    }

    /// An upper hull with one point per x, strictly concave, covers(p) means p lies on or below it
    struct Chain{
        map<int64_t, int64_t> pts;

        static Point at(map<int64_t, int64_t>::const_iterator it){
            return Point(it->first, it->second);
        }

        bool covers(const Point& p) const {
            auto it = pts.lower_bound(p.x);
            if (it == pts.end()) return false;
            if (it->first == p.x) return p.y <= it->second;
            if (it == pts.begin()) return false;
            return cross(at(prev(it)), at(it), p) <= 0;
        }

        void add(const Point& p){
            if (covers(p)) return;
            pts[p.x] = p.y;
            auto it = pts.find(p.x);

            while (next(it) != pts.end() && next(next(it)) != pts.end() && cross(p, at(next(it)), at(next(next(it)))) >= 0){
                pts.erase(next(it));
            }
            while (it != pts.begin() && prev(it) != pts.begin() && cross(at(prev(prev(it))), at(prev(it)), p) >= 0){
                pts.erase(prev(it));
            }
        }
    };

    Chain upper, lower;  /// lower holds the points mirrored by y -> -y
};

void solve(vector<Point> polygon, bool first){
    /// The polygon repeats its first vertex at the end
    polygon.pop_back();
    vector<Point> hull = get_convex_hull(polygon);
    vector<Point> sorted_hull = hull;
    sort(sorted_hull.begin(), sorted_hull.end());
    int i = 0;
    while (!binary_search(sorted_hull.begin(), sorted_hull.end(), polygon[i])) i++;
    int start = find(hull.begin(), hull.end(), polygon[i]) - hull.begin();

    if (!first) puts("");
    int n = hull.size();
    for (int i = 0; i <= n; i++){
        const Point& p = hull[(start + i) % n];
        printf("%lld, %lld\n", (long long)p.x, (long long)p.y);
    }
}

int main(){
    char line[256];
    bool first = true;
    vector<Point> polygon;
    while (fgets(line, sizeof(line), stdin)){
        long long x, y;
        if (sscanf(line, "%lld , %lld", &x, &y) == 2) polygon.push_back(Point(x, y));
        else if (!polygon.empty()){
            solve(polygon, first);
            first = false, polygon.clear();
        }
    }
    if (!polygon.empty()) solve(polygon, first);
    return 0;
}
