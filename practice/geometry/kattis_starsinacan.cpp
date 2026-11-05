/***
 *
 * 3D Convex Hull
 * Triangulated convex hull of integer points in 3D, with its volume and surface area
 *
 * Complexity: expected O(n log n) time, memory peaks near 150 MB for 1e5 points that are all vertices
 *
 * ConvexHull3D hull(points): the hull of points, faces empty when every point lies on one plane
 * hull.area(): surface area
 * hull.volume(): volume
 * hull.faces: triangles {a, b, c} as indices into points, counter-clockwise seen from outside
 *
 * Notes:
 *   - |x|, |y|, |z| <= 1e9, every orientation test is exact in __int128 (64-bit GCC or Clang)
 *   - No input point lies strictly outside the plane of any face, and no face has zero area
 *   - Duplicate points are allowed, a vertex is any one of its copies
 *   - A flat facet comes as coplanar triangles, maybe with vertices inside it or on its edges; count facets by plane
 *   - The insertion order is random, so the triangulation of a flat facet can change between runs
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point3D{
    long long x, y, z;
};

struct ConvexHull3D{
    struct Face{
        array<int, 3> v, next;  /// next[i] is the face across the edge v[i] -> v[(i + 1) % 3]
        long long nx, ny, nz;
        vector<int> outside;  /// unadded points strictly above this face, emptied when the face dies
        bool dead = false;
    };

    vector<Point3D> points;
    vector<array<int, 3>> faces;

    ConvexHull3D(const vector<Point3D>& points) : points(points){
        int n = points.size(), a = 0, b = 1, c = 0, d = 0;
        auto same = [&](int i){
            return points[i].x == points[a].x && points[i].y == points[a].y && points[i].z == points[a].z;
        };
        auto collinear = [&](int i){
            Face f = make_face(a, b, i);
            return !f.nx && !f.ny && !f.nz;
        };
        while (b < n && same(b)) b++;
        while (b < n && c < n && collinear(c)) c++;
        if (b >= n || c >= n) return;

        Face abc = make_face(a, b, c);
        while (d < n && side(abc, d) == 0) d++;
        if (d == n) return;

        if (side(abc, d) > 0) swap(b, c);
        vector<Face> all;
        int base[4][3] = {{a, b, c}, {a, d, b}, {b, d, c}, {c, d, a}};
        int across[4][3] = {{1, 2, 3}, {3, 2, 0}, {1, 3, 0}, {2, 1, 0}};
        for (int f = 0; f < 4; f++){
            all.push_back(make_face(base[f][0], base[f][1], base[f][2]));
            all[f].next = {across[f][0], across[f][1], across[f][2]};
        }

        vector<vector<int>> sees(n);
        vector<int> order;
        for (int p = 0; p < n; p++){
            if (p == a || p == b || p == c || p == d) continue;
            order.push_back(p);
            for (int f = 0; f < 4; f++){
                if (side(all[f], p) > 0) all[f].outside.push_back(p), sees[p].push_back(f);
            }
        }
        /// Random insertion order is what makes the expected total conflict work O(n log n)
        static mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
        shuffle(order.begin(), order.end(), rng);

        vector<int> start(n), stamp(n, -1), visible, fresh;
        for (int p : order){
            visible.clear();
            for (int f : sees[p]){
                if (!all[f].dead) all[f].dead = true, visible.push_back(f);
            }
            sees[p] = vector<int>();
            if (visible.empty()) continue;

            fresh.clear();
            for (int f : visible){
                for (int i = 0; i < 3; i++){
                    int g = all[f].next[i];
                    if (all[g].dead) continue;

                    int u = all[f].v[i], w = all[f].v[(i + 1) % 3], h = all.size();
                    all.push_back(make_face(u, w, p));
                    all[h].next[0] = g;
                    for (int j = 0; j < 3; j++){
                        if (all[g].v[j] == w) all[g].next[j] = h;
                    }
                    start[u] = h;
                    fresh.push_back(h);

                    for (const vector<int>* list : {&all[f].outside, &all[g].outside}){
                        for (int q : *list){
                            if (stamp[q] == h || side(all[h], q) <= 0) continue;

                            stamp[q] = h;
                            all[h].outside.push_back(q);
                            sees[q].push_back(h);
                        }
                    }
                }
            }

            for (int h : fresh) all[h].next[1] = start[all[h].v[1]];
            for (int h : fresh) all[all[h].next[1]].next[2] = h;
            for (int f : visible) all[f].outside = vector<int>();
        }

        for (const Face& f : all){
            if (!f.dead) faces.push_back(f.v);
        }
    }

    long double area() const{
        long double res = 0;
        for (auto& f : faces){
            Face g = make_face(f[0], f[1], f[2]);
            res += sqrtl((long double)g.nx * g.nx + (long double)g.ny * g.ny + (long double)g.nz * g.nz);
        }

        return res / 2;
    }

    long double volume() const{
        __int128 res = 0;
        for (auto& f : faces) res -= side(make_face(f[0], f[1], f[2]), faces[0][0]);

        return (long double)res / 6;
    }

    Face make_face(int a, int b, int c) const{
        const Point3D &p = points[a], &q = points[b], &r = points[c];
        long long ux = q.x - p.x, uy = q.y - p.y, uz = q.z - p.z, vx = r.x - p.x, vy = r.y - p.y, vz = r.z - p.z;
        Face f;
        f.v = {a, b, c};
        f.nx = uy * vz - uz * vy, f.ny = uz * vx - ux * vz, f.nz = ux * vy - uy * vx;
        return f;
    }

    /// Six times the signed volume of the face's triangle and q, positive when q is outside
    __int128 side(const Face& f, int q) const{
        const Point3D &p = points[f.v[0]], &r = points[q];
        return (__int128)f.nx * (r.x - p.x) + (__int128)f.ny * (r.y - p.y) + (__int128)f.nz * (r.z - p.z);
    }
};

/***
 *
 * Circle Geometry
 * Intersections, tangents, enclosing circles, covering and areas of circles in long double
 *
 * Complexity:
 *   - intersections, tangents, circumcircle, circle_circle_area: O(1)
 *   - circle_polygon_area: O(n)
 *   - minimum_enclosing_circle: expected O(n)
 *   - circle_union_area, max_circle_cover: O(n^2 log n)
 *
 * Circle{center, r}: a circle with r >= 0
 * circle_circle_area(a, b): area of the intersection of two circles
 * circle_circle_intersection(a, b): 0, 1 (tangent) or 2 points, none for concentric circles (identical ones included)
 * circle_line_intersection(c, a, b): where the line through a != b meets c, 0, 1 or 2 points in the a -> b direction
 * circle_polygon_area(c, poly): area of c intersected with a simple polygon in either orientation
 * circle_union_area(circles): area of the union
 * circumcircle(a, b, c): the circle through three non-collinear points
 * common_tangents(a, b): every common tangent as its tangency points {on a, on b}, r > 0 for both, outer ones first
 * max_circle_cover(points, r): {most points a radius r circle covers, a center achieving it}, boundary points count
 * minimum_enclosing_circle(points): smallest circle containing every point, points non-empty
 * tangents_from_point(c, p): tangency points of the lines through p touching c (r > 0), none if p is inside, p if on c
 *
 * Notes:
 *   - Tangency uses a tolerance EPS relative to r, so a circle with r = 0 meets a line only through its center exactly
 *   - Of 2 intersection or tangency points, the first lies left of the line from the first center to the second or p
 *   - common_tangents gives 0 lines nested, 1 internally tangent, 2 overlapping, 3 externally tangent, 4 apart
 *   - common_tangents gives none for identical circles, a line whose points coincide is perpendicular to the centers
 *   - Verified with |coordinates|, radii <= 1e4
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

using Float = long double;

const Float EPS = 1e-12L;
const Float PI = acosl(-1.0L);

struct Point{
    Float x, y;

    Point operator+(const Point& p) const{
        return {x + p.x, y + p.y};
    }

    Point operator-(const Point& p) const{
        return {x - p.x, y - p.y};
    }

    Point operator*(Float k) const{
        return {x * k, y * k};
    }

    Point operator/(Float k) const{
        return {x / k, y / k};
    }

    Float dot(const Point& p) const{
        return x * p.x + y * p.y;
    }

    Float cross(const Point& p) const{
        return x * p.y - y * p.x;
    }

    Float norm2() const{
        return x * x + y * y;
    }

    Float norm() const{
        return sqrtl(norm2());
    }

    Point perp() const{
        return {-y, x};
    }
};

struct Circle{
    Point center;
    Float r;
};

Float circle_circle_area(const Circle& a, const Circle& b){
    Float d = (a.center - b.center).norm(), r = min(a.r, b.r);
    if (d >= a.r + b.r) return 0;
    if (d <= fabsl(a.r - b.r)) return PI * r * r;

    Float alpha = acosl(clamp((d * d + a.r * a.r - b.r * b.r) / (2 * d * a.r), (Float)-1, (Float)1));
    Float beta = acosl(clamp((d * d + b.r * b.r - a.r * a.r) / (2 * d * b.r), (Float)-1, (Float)1));
    return a.r * a.r * (alpha - sinl(2 * alpha) / 2) + b.r * b.r * (beta - sinl(2 * beta) / 2);
}

vector<Point> circle_circle_intersection(const Circle& a, const Circle& b){
    Point d = b.center - a.center;
    Float d2 = d.norm2(), sum = a.r + b.r, dif = a.r - b.r, tol = EPS * sum * sum;

    if (d2 == 0 || d2 > sum * sum + tol || d2 < dif * dif - tol) return {};

    Float p = (d2 + a.r * a.r - b.r * b.r) / (2 * d2);
    Point mid = a.center + d * p;
    if (fabsl(d2 - sum * sum) <= tol || fabsl(d2 - dif * dif) <= tol) return {mid};

    Point h = d.perp() * sqrtl(max(a.r * a.r - p * p * d2, (Float)0) / d2);
    return {mid + h, mid - h};
}

vector<Point> circle_line_intersection(const Circle& c, const Point& a, const Point& b){
    Point ab = b - a;
    Float len2 = ab.norm2(), s = ab.cross(c.center - a);
    Float h2 = c.r * c.r * len2 - s * s, tol = EPS * c.r * c.r * len2;

    if (h2 < -tol) return {};

    Point foot = a + ab * (ab.dot(c.center - a) / len2);
    if (h2 <= tol) return {foot};

    Point h = ab * (sqrtl(h2) / len2);
    return {foot - h, foot + h};
}

Float circle_polygon_area(const Circle& c, const vector<Point>& poly){
    Float half_r2 = c.r * c.r / 2, res = 0;
    auto angle = [](const Point& p, const Point& q){ return atan2l(p.cross(q), p.dot(q)); };

    /// Signed area of the circle intersected with the triangle (center, p, q), center at the origin
    auto sector = [&](const Point& p, const Point& q) -> Float{
        Point d = q - p;
        Float len2 = d.norm2();
        if (len2 == 0) return 0;

        Float a = d.dot(p) / len2, b = (p.norm2() - c.r * c.r) / len2, det = a * a - b;
        if (det <= 0) return angle(p, q) * half_r2;

        Float s = max((Float)0, -a - sqrtl(det)), t = min((Float)1, -a + sqrtl(det));
        if (t < 0 || s >= 1) return angle(p, q) * half_r2;

        Point u = p + d * s, v = p + d * t;
        return angle(p, u) * half_r2 + u.cross(v) / 2 + angle(v, q) * half_r2;
    };

    for (int i = 0, n = poly.size(); i < n; i++) res += sector(poly[i] - c.center, poly[(i + 1) % n] - c.center);
    return fabsl(res);
}

Float circle_union_area(const vector<Circle>& circles){
    int n = circles.size();
    Float res = 0;

    for (int i = 0; i < n; i++){
        const Circle& c = circles[i];
        bool hidden = c.r <= 0;
        vector<pair<Float, Float>> covered;
        for (int j = 0; j < n && !hidden; j++){
            const Circle& o = circles[j];
            Point v = o.center - c.center;
            Float d = v.norm();

            /// Of identical circles only the first one counts, otherwise each would hide the other
            if (j != i && d <= o.r - c.r && (o.r > c.r || j < i)) hidden = true;
            if (j == i || d >= c.r + o.r || d <= fabsl(c.r - o.r)) continue;

            Float mid = atan2l(v.y, v.x), half = acosl(clamp((c.r * c.r + d * d - o.r * o.r) / (2 * c.r * d), (Float)-1, (Float)1));
            Float lo = mid - half, hi = mid + half;
            if (lo < 0) lo += 2 * PI, hi += 2 * PI;
            if (hi <= 2 * PI) covered.push_back({lo, hi});
            else covered.push_back({lo, 2 * PI}), covered.push_back({0, hi - 2 * PI});
        }
        if (hidden) continue;

        /// Green's theorem: the arc from angle s to t adds half the integral of x dy - y dx along it
        auto arc = [&](Float s, Float t){
            return (c.r * c.r * (t - s) + c.center.x * c.r * (sinl(t) - sinl(s)) - c.center.y * c.r * (cosl(t) - cosl(s))) / 2;
        };

        sort(covered.begin(), covered.end());
        Float from = 0;
        for (auto [lo, hi] : covered){
            if (lo > from) res += arc(from, lo);
            from = max(from, hi);
        }
        if (from < 2 * PI) res += arc(from, 2 * PI);
    }

    return res;
}

Circle circumcircle(const Point& a, const Point& b, const Point& c){
    Point u = b - a, v = c - a;
    Point center = a + (v * u.norm2() - u * v.norm2()).perp() / (v.cross(u) * 2);
    return {center, max({(center - a).norm(), (center - b).norm(), (center - c).norm()})};
}

vector<pair<Point, Point>> common_tangents(const Circle& a, const Circle& b){
    vector<pair<Point, Point>> res;
    Point d = b.center - a.center;
    Float d2 = d.norm2();

    if (d2 == 0) return res;

    /// A negated radius puts the tangency point on the far side of b, turning outer tangents into inner ones
    for (Float rb : {b.r, -b.r}){
        Float dr = a.r - rb, h2 = d2 - dr * dr;
        if (h2 < -EPS * d2) continue;

        Float h = h2 <= EPS * d2 ? 0 : sqrtl(h2);
        for (int sign : {1, -1}){
            Point v = (d * dr + d.perp() * (h * sign)) / d2;
            res.push_back({a.center + v * a.r, b.center + v * rb});
            if (h == 0) break;
        }
    }

    return res;
}

pair<int, Point> max_circle_cover(const vector<Point>& points, Float r){
    int n = points.size(), best = 0;
    Point center{0, 0};

    for (int i = 0; i < n; i++){
        int count = 0;
        vector<pair<Float, int>> events;
        for (int j = 0; j < n; j++){
            Point v = points[j] - points[i];
            Float d2 = v.norm2();
            if (d2 == 0){
                count++;
                continue;
            }
            if (d2 > 4 * r * r * (1 + EPS)) continue;

            /// Centers at distance r from points[i] that also cover points[j], widened by EPS for boundary ties
            Float mid = atan2l(v.y, v.x), half = atan2l(sqrtl(max((Float)0, 4 * r * r - d2)), sqrtl(d2));
            Float lo = mid - half - EPS, hi = mid + half + EPS;
            if (lo < 0) lo += 2 * PI, hi += 2 * PI;
            events.push_back({lo, -1});
            events.push_back({min(hi, 2 * PI), 1});
            if (hi > 2 * PI) events.push_back({0, -1}), events.push_back({hi - 2 * PI, 1});
        }

        /// -1 opens an interval, so at equal angles it sorts before the closing +1
        sort(events.begin(), events.end());
        Float at = 0;
        int most = count;
        for (auto [angle, type] : events){
            count -= type;
            if (count > most) most = count, at = angle;
        }
        if (most > best) best = most, center = points[i] + Point{cosl(at), sinl(at)} * r;
    }

    return {best, center};
}

Circle minimum_enclosing_circle(vector<Point> points){
    assert(!points.empty());
    shuffle(points.begin(), points.end(), mt19937(chrono::steady_clock::now().time_since_epoch().count()));

    int n = points.size();
    Circle c{points[0], 0};
    auto outside = [&](const Point& p){ return (p - c.center).norm() > c.r * (1 + EPS); };
    for (int i = 1; i < n; i++){
        if (!outside(points[i])) continue;
        c = {points[i], 0};
        for (int j = 0; j < i; j++){
            if (!outside(points[j])) continue;
            /// A rounded midpoint radius could leave pi outside and pass a degenerate triple to circumcircle (NaN)
            Point mid = (points[i] + points[j]) / 2;
            c = {mid, max((mid - points[i]).norm(), (mid - points[j]).norm())};
            for (int k = 0; k < j; k++){
                if (outside(points[k])) c = circumcircle(points[i], points[j], points[k]);
            }
        }
    }

    return c;
}

vector<Point> tangents_from_point(const Circle& c, const Point& p){
    Point d = p - c.center;
    Float d2 = d.norm2(), r2 = c.r * c.r;
    if (d2 < r2 * (1 - EPS)) return {};
    if (d2 <= r2 * (1 + EPS)) return {p};

    Point foot = c.center + d * (r2 / d2), h = d.perp() * (c.r * sqrtl(d2 - r2) / d2);
    return {foot + h, foot - h};
}

/***
 *
 * Fast I/O
 * Buffered stdin and stdout, about 5x faster than scanf for integers, interactive problems included
 *
 * Complexity: O(1) amortized per character, one system call per buffer of up to 8192 bytes
 *
 * fio::flush(): writes the pending output now, which also happens at normal exit
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ..., returns the count read before the input ended
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - Reading into a char array does not compile: read a string, or a grid as vector<string>
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
 *   - Interactive problems need no flush(): pending output is written each time the input buffer runs empty
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>
#include <unistd.h>  /// read() and write(), which MinGW's unistd.h also declares, so Windows judges compile them too

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// write, not fwrite: fwrite strands bytes in stdout's buffer, and fflush after it costs a second call, 7% slower
    void flush(){
        int done = 0;
        while (done < outptr){
            int len = ::write(1, outbuf + done, outptr - done);
            if (len <= 0) break;
            done += len;
        }
        outptr = 0;
    }

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            /// An interactive judge answers only after it has read our query, so the query leaves before we wait
            flush();

            /// read returns the bytes available, fread waits for all 8192 or EOF: a deadlock on an interactive pipe
            inptr = 0, buf_len = ::read(0, inbuf, BUF_SIZE);
            if (buf_len <= 0) return EOF;
        }

        return (unsigned char)inbuf[inptr++];
    }

    /// Catches every integer the overloads below miss, so __int128 needs no is_integral specialization (breaks C++20)
    template <typename T>
    bool read_one(T &x){
        int c = ' ', neg = 0;

        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// A local, since an int x may alias inptr or buf_len and writing it per digit forces their reload from memory
        T v = 0;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (; isdigit(c); c = read_char()){
            v = v * 10 - (c - '0');
        }
        x = neg ? v : -v;

        return true;
    }

    bool read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        if (!read_token()) return false;
        x = strtof(tmpbuf, nullptr);
        return true;
    }

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
        const unsigned char *s = (const unsigned char *)tmpbuf + (*tmpbuf == '-' || *tmpbuf == '+');
        unsigned long long m = 0;
        int digits = 0, e = 0;
        for (; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++;
        if (*s == '.'){
            for (s++; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++, e--;
        }

        /// No exponent digits after the e equals strtod, which then ignores the e
        if (*s == 'e' || *s == 'E'){
            int sign = s[1] == '-' ? -1 : 1, ex = 0;
            s += 1 + (s[1] == '-' || s[1] == '+');
            for (; isdigit(*s) && ex < 1000; s++) ex = ex * 10 + (*s - '0');
            e += sign * ex;
        }

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
        return true;
    }

    bool read_one(long double &x){
        if (!read_token()) return false;
        x = strtold(tmpbuf, nullptr);
        return true;
    }

    bool read_one(string &s){
        int c = ' ';

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !isspace(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T> &v){
        for (auto &&x: v){
            if (!read_one(x)) return false;
        }

        return true;
    }

    template <typename T>
    int read(T &x){
        return read_one(x);
    }

    template <typename T, typename ...Args>
    int read(T &x, Args &...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
        return true;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(char c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    /// Integers, like read_one: the plain overloads below win for every other type
    template <typename T>
    void write_one(T x){
        int l = 0;

        if (x < 0){
            /// Peel one digit while still negative so that -x cannot overflow for the minimum value
            write_char('-');
            tmpbuf[l++] = '0' - x % 10;
            x = -(x / 10);
        }

        while (x || !l){
            tmpbuf[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmpbuf[--l]);
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char *s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    /// Without it a char * picks the integer template, an exact match that beats the const char * overload
    void write_one(char *s){
        write_one((const char *)s);
    }

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
        static const long long p10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
            10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 1000000000000000,
            10000000000000000, 100000000000000000, 1000000000000000000};

        if (f.digits >= 0 && f.digits <= 18){
            long double scaled = fabsl(f.value) * p10[f.digits];

            /// Rounding never moves the product across k + 0.5, exact below 2^52, so only an exact k + 0.5 needs printf
            if (scaled < 4503599627370496.0L && scaled - (long long)scaled != 0.5L){
                long long q = scaled + 0.5L;
                if (signbit(f.value)) write_char('-');
                write_one(q / p10[f.digits]);
                if (!f.digits) return;

                long long r = q % p10[f.digits];
                for (int j = f.digits - 1; j >= 0; j--) tmpbuf[j] = '0' + r % 10, r /= 10;
                write_char('.');
                for (int j = 0; j < f.digits; j++) write_char(tmpbuf[j]);
                return;
            }
        }

        int len = snprintf(tmpbuf, sizeof(tmpbuf), "%.*Lf", f.digits, f.value);

        if (len < (int)sizeof(tmpbuf)){
            for (int j = 0; j < len; j++) write_char(tmpbuf[j]);
            return;
        }

        string s(len + 1, '\0');
        snprintf(&s[0], len + 1, "%.*Lf", f.digits, f.value);
        s.resize(len);
        write_one(s);
    }

    template <typename T>
    void write_one(const vector<T> &v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    template <typename T>
    char separator_after(const T &){
        return ' ';
    }

    template <typename T>
    char separator_after(const vector<T> &){
        return '\n';
    }

    template <typename T>
    void write(const T &x){
        write_one(x);
        write_char('\n');
    }

    template <typename T, typename ...Args>
    void write(const T &x, const Args &...args){
        write_one(x);
        write_char(separator_after(x));
        write(args...);
    }

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }
}

using namespace fio;

int n;
vector<Point3D> stars;

Float volume(const ConvexHull3D& hull, const array<int, 3>& f){
    int i;
    Float h = 0, norm, len, nx, ny, nz, ux, uy, uz, vx, vy, vz;

    auto g = hull.make_face(f[0], f[1], f[2]);
    norm = sqrtl((Float)g.nx * g.nx + (Float)g.ny * g.ny + (Float)g.nz * g.nz);
    nx = g.nx / norm, ny = g.ny / norm, nz = g.nz / norm;

    if (fabsl(nx) <= fabsl(ny) && fabsl(nx) <= fabsl(nz)) ux = 0, uy = nz, uz = -ny;
    else if (fabsl(ny) <= fabsl(nz)) ux = -nz, uy = 0, uz = nx;
    else ux = ny, uy = -nx, uz = 0;
    len = sqrtl(ux * ux + uy * uy + uz * uz);
    ux /= len, uy /= len, uz /= len;
    vx = ny * uz - nz * uy, vy = nz * ux - nx * uz, vz = nx * uy - ny * ux;

    vector<Point> points(n);
    for (i = 0; i < n; i++){
        h = max(h, -hull.side(g, i) / norm);
        points[i] = {ux * stars[i].x + uy * stars[i].y + uz * stars[i].z, vx * stars[i].x + vy * stars[i].y + vz * stars[i].z};
    }

    Circle c = minimum_enclosing_circle(points);
    return PI * c.r * c.r * h;
}

int main(){
    int i;
    Float res = 1e30;

    read(n);
    stars.resize(n);
    for (i = 0; i < n; i++) read(stars[i].x, stars[i].y, stars[i].z);

    ConvexHull3D hull(stars);
    for (auto& f: hull.faces) res = min(res, volume(hull, f));

    write(fixed(res, 8));
    return 0;
}
