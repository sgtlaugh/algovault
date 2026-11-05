/***
 *
 * Convex Hull
 * Monotone chain hull, rotating calipers, Minkowski sum, O(log n) queries on a hull and an online insertion hull
 *
 * Complexity:
 *   - get_convex_hull: O(n log n)
 *   - farthest_pair, is_convex, width: O(n)
 *   - minkowski_sum: O(n + m)
 *   - extreme_vertex, line_hull, tangents: O(log n)
 *   - DynamicHull: add amortized O(log n), contains O(log n), hull O(h)
 *
 * extreme_vertex(hull, dir): index of a vertex maximizing dot(vertex, dir), n >= 1, |dir.x|, |dir.y| <= 2e9
 * farthest_pair(hull): two vertices at maximum distance, dist2 of them is the exact squared diameter, n >= 1
 * get_convex_hull(points): the strict hull
 * is_convex(polygon): whether polygon is convex
 * line_hull(a, b, hull): how the line through a != b meets a hull with n >= 3
 * minkowski_sum(a, b): the strict hull of {p + q}, a and b convex counter-clockwise polygons
 * tangents(hull, q): {i, j}, the hull lies left of the ray q -> hull[i] and right of q -> hull[j], q strictly outside
 * width(hull): minimum distance between two parallel lines enclosing the hull, 0 when n <= 2
 * DynamicHull dh: an empty hull taking points online
 * dh.add(p): inserts p
 * dh.contains(p): whether p lies inside the hull or on its boundary
 * dh.hull(): the current hull as get_convex_hull returns it
 *
 * Notes:
 *   - |x|, |y| <= 1e9, so every cross product, dot product and squared distance fits in int64_t
 *   - Input with k fixed decimal places: take llround(x * 10^k), exact, where a cast truncates 0.29 * 100 to 28
 *   - A hull argument is a strictly convex polygon in counter-clockwise order, as get_convex_hull returns
 *   - line_hull gives {-1, -1} for no intersection, {i, -1} touching only vertex i, {i, i} containing edge (i, i + 1)
 *   - Else line_hull gives {i, j}, crossing edge (i, i + 1) then (j, j + 1) going a -> b, a vertex k as edge (k, k + 1)
 *   - minkowski_sum takes any starting vertex and collinear vertices, its coordinates reach 2e9
 *   - tangents returns either endpoint when q is collinear with an edge
 *   - get_convex_hull becomes O(n) with a radix sort
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    int64_t x, y;

    Point() {}
    Point(int64_t x, int64_t y) : x(x), y(y) {}

    inline bool operator < (const Point &p) const{
        return ((x < p.x) || (x == p.x && y < p.y));
    }

    bool operator == (const Point& p) const{
        return x == p.x && y == p.y;
    }

    Point operator + (const Point& p) const{
        return Point(x + p.x, y + p.y);
    }

    Point operator - (const Point& p) const{
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

    /// Same-sign turns still allow stars, a convex polygon reverses its x and its y direction exactly twice each
    return sign != 0 && flips[0] == 2 && flips[1] == 2;
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

/// Seen from an outside point the hull spans less than 180 degrees, so the angular order behaves like a projection
pair<int, int> tangents(const vector<Point>& hull, const Point& q){
    int n = hull.size();
    int right = extreme_search(n, [&](int i, int j){ return sgn(cross(q, hull[i], hull[j])); });
    int left = extreme_search(n, [&](int i, int j){ return sgn(cross(q, hull[j], hull[i])); });
    return {right, left};
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

struct DynamicHull{
    /// An upper hull with one point per x, strictly concave, covers(p) means p lies on or below it
    struct Chain{
        map<int64_t, int64_t> pts;

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

        bool covers(const Point& p) const{
            auto it = pts.lower_bound(p.x);
            if (it == pts.end()) return false;
            if (it->first == p.x) return p.y <= it->second;
            if (it == pts.begin()) return false;
            return cross(at(prev(it)), at(it), p) <= 0;
        }

        static Point at(map<int64_t, int64_t>::const_iterator it){
            return Point(it->first, it->second);
        }
    };

    Chain upper, lower;  /// lower holds the points mirrored by y -> -y

    void add(const Point& p){
        upper.add(p);
        lower.add(Point(p.x, -p.y));
    }

    bool contains(const Point& p) const{
        return upper.covers(p) && lower.covers(Point(p.x, -p.y));
    }

    vector<Point> hull() const{
        vector<Point> res;

        for (auto& [x, y] : lower.pts) res.push_back(Point(x, -y));

        for (auto it = upper.pts.rbegin(); it != upper.pts.rend(); it++){
            Point p(it->first, it->second);
            if (!(p == res.back()) && !(p == res[0])) res.push_back(p);
        }

        return res;
    }
};

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

int main(){
    int n, i, m;
    long double r, d;

    read(n, r);
    vector<Point> points(n);
    for (i = 0; i < n; i++) read(points[i].x, points[i].y);

    vector<Point> hull = get_convex_hull(points);
    m = hull.size(), d = r;
    Point o(0, 0);
    for (i = 0; i < m && m >= 3; i++){
        Point u = hull[i], v = hull[(i + 1) % m];
        d = min(d, max((long double)0, cross(u, v, o) / sqrtl(dist2(u, v))));
    }
    if (m < 3) d = 0;

    long double res = (r * r * acosl(d / r)) - (d * sqrtl((r * r) - (d * d)));
    write(fixed(res, 12));

    return 0;
}
