/***
 *
 * Half-Plane Intersection
 * Intersects half-planes, each the left side of a directed line, inside the square |x|, |y| <= box
 *
 * Complexity: O(n log n)
 *
 * Halfplane(a, b): the points p with cross(b - a, p - a) >= 0, boundary included, a != b
 * IntHalfplane(ax, ay, bx, by): the same half-plane for integer a = (ax, ay) != b = (bx, by), |coordinates| <= 1e9
 * halfplane_intersection(h, box): vertices of the intersection inside the box, counter-clockwise, strictly convex
 *
 * Notes:
 *   - h is all Halfplane or all IntHalfplane; with IntHalfplane, box is an integer in [1, 1e9]
 *   - The starting vertex is unspecified and no vertex repeats
 *   - Empty if the intersection is empty or has zero area (a segment or a single point)
 *   - A bounded intersection comes out exact when box exceeds all its vertex coordinates
 *   - IntHalfplane decides exactly in __int128, rounds only the output, and is 3x faster: use it for integer input
 *   - Halfplane is long double with an absolute EPS = 1e-9, a vertex that close to a boundary is on it, so box <= 1e9
 *   - Halfplane merges directions within the double rounding of the points, so a line given twice by decimals merges
 *   - Halfplane never merges distinct integer directions with |x|, |y| <= 1e6
 *   - Verified against brute force: Halfplane for integer |x|, |y| <= 1e4 and 0.1-step decimals, IntHalfplane exactly
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long double EPS = 1e-9;

struct Point{
    long double x, y;

    explicit Point(long double x = 0, long double y = 0) : x(x), y(y) {}

    Point operator+(const Point& p) const{
        return Point(x + p.x, y + p.y);
    }

    Point operator-(const Point& p) const{
        return Point(x - p.x, y - p.y);
    }

    Point operator*(long double k) const{
        return Point(x * k, y * k);
    }
};

long double cross(const Point& a, const Point& b){
    return a.x * b.y - a.y * b.x;
}

struct Halfplane{
    Point p, pq;
    long double len;

    Halfplane(const Point& a, const Point& b) : p(a), pq(b - a), len(hypotl(pq.x, pq.y)) {}

    bool out(const Point& r) const{
        return cross(pq, r - p) < -EPS * len;
    }

    /// Angle in [0, pi), tested on the exact direction so -0.0 and parallel lines stay consistent
    bool upper() const{
        return pq.y > 0 || (pq.y == 0 && pq.x > 0);
    }
};

struct IntHalfplane{
    static constexpr long long LIMIT = 1000000000;

    long long px, py, dx, dy;

    IntHalfplane(long long ax, long long ay, long long bx, long long by) : px(ax), py(ay), dx(bx - ax), dy(by - ay){
        assert((dx != 0 || dy != 0) && max({llabs(ax), llabs(ay), llabs(bx), llabs(by)}) <= LIMIT);
    }

    long long cross_dir(const IntHalfplane& h) const{
        return dx * h.dy - dy * h.dx;
    }

    long long side(long long x, long long y) const{
        return dx * (y - py) - dy * (x - px);
    }

    /// Needs a.cross_dir(b) > 0, or the sign flips; with |coordinates| <= LIMIT both products stay below 2^126
    bool out(const IntHalfplane& a, const IntHalfplane& b) const{
        long long det = a.cross_dir(b), t = b.side(a.px, a.py);
        return (__int128)side(a.px, a.py) * det + (__int128)cross_dir(a) * t <= 0;
    }

    bool upper() const{
        return dy > 0 || (dy == 0 && dx > 0);
    }
};

Point line_intersection(const Halfplane& s, const Halfplane& t){
    long double alpha = cross(t.p - s.p, t.pq) / cross(s.pq, t.pq);
    return s.p + s.pq * alpha;
}

vector<Point> halfplane_intersection(vector<Halfplane> h, long double box){
    Point corners[4] = {Point(box, -box), Point(box, box), Point(-box, box), Point(-box, -box)};

    for (int i = 0; i < 4; i++) h.emplace_back(corners[i], corners[(i + 1) % 4]);

    sort(h.begin(), h.end(), [](const Halfplane& a, const Halfplane& b){
        if (a.upper() != b.upper()) return a.upper();
        return cross(a.pq, b.pq) > 0;
    });

    /// Of lines with directions equal up to double rounding of their points, only the innermost constrains anything
    auto size = [](const Halfplane& a){ return fabsl(a.p.x) + fabsl(a.p.y) + fabsl(a.pq.x) + fabsl(a.pq.y); };
    /// The dot product, not upper(), splits antiparallel lines, so twins rounded across angle pi still merge
    auto same_direction = [&](const Halfplane& a, const Halfplane& b){
        return a.pq.x * b.pq.x + a.pq.y * b.pq.y > 0 && fabsl(cross(a.pq, b.pq)) <= 8 * DBL_EPSILON * (size(a) * b.len + size(b) * a.len);
    };
    auto merge = [](Halfplane& kept, const Halfplane& hp){
        if (cross(kept.pq, hp.p - kept.p) > 0) kept = hp;
    };
    vector<Halfplane> lines;
    for (auto& hp : h){
        if (!lines.empty() && same_direction(lines.back(), hp)) merge(lines.back(), hp);
        else lines.push_back(hp);
    }

    /// A direction rounded to just below angle 0 sorts last, next to its twin at the front
    if (same_direction(lines.back(), lines[0])){
        merge(lines[0], lines.back());
        lines.pop_back();
    }

    deque<Halfplane> dq;
    for (auto& hp : lines){
        while (dq.size() > 1 && hp.out(line_intersection(dq[dq.size() - 1], dq[dq.size() - 2]))) dq.pop_back();
        while (dq.size() > 1 && hp.out(line_intersection(dq[0], dq[1]))) dq.pop_front();

        /// Everything between two opposite directions was cut away, so nothing is left with positive area
        if (!dq.empty() && cross(dq.back().pq, hp.pq) <= 0) return {};
        dq.push_back(hp);
    }

    while (dq.size() > 2 && dq[0].out(line_intersection(dq[dq.size() - 1], dq[dq.size() - 2]))) dq.pop_back();
    while (dq.size() > 2 && dq.back().out(line_intersection(dq[0], dq[1]))) dq.pop_front();
    if (dq.size() < 3) return {};

    /// Lines through a common vertex leave zero-length edges, a degenerate intersection collapses to under 3 vertices
    auto same = [](const Point& a, const Point& b){ return fabsl(a.x - b.x) <= EPS && fabsl(a.y - b.y) <= EPS; };
    vector<Point> poly;
    for (int i = 0, k = dq.size(); i < k; i++){
        Point v = line_intersection(dq[i], dq[(i + 1) % k]);
        if (poly.empty() || !same(poly.back(), v)) poly.push_back(v);
    }
    while (poly.size() > 1 && same(poly.back(), poly[0])) poly.pop_back();

    if (poly.size() < 3) return {};
    return poly;
}

/// Exact numerators leave two roundings per coordinate, the conversion and the division, so about 1 ulp
Point line_intersection(const IntHalfplane& s, const IntHalfplane& t){
    long long det = s.cross_dir(t), alpha = t.side(s.px, s.py);
    __int128 x = (__int128)s.px * det + (__int128)s.dx * alpha, y = (__int128)s.py * det + (__int128)s.dy * alpha;
    return Point((long double)x / det, (long double)y / det);
}

vector<Point> halfplane_intersection(vector<IntHalfplane> h, long long box){
    assert(1 <= box && box <= IntHalfplane::LIMIT);
    long long cx[4] = {box, box, -box, -box}, cy[4] = {-box, box, box, -box};
    for (int i = 0; i < 4; i++) h.emplace_back(cx[i], cy[i], cx[(i + 1) % 4], cy[(i + 1) % 4]);

    /// Splitting by upper() first, then sorting each half by cross, is 13% faster than one comparator at n = 5e4
    auto mid = partition(h.begin(), h.end(), [](const IntHalfplane& a){ return a.upper(); });
    auto by_angle = [](const IntHalfplane& a, const IntHalfplane& b){ return a.cross_dir(b) > 0; };
    sort(h.begin(), mid, by_angle);
    sort(mid, h.end(), by_angle);

    /// Antiparallel lines are never neighbors, as the box puts its directions pi / 2 and pi between them
    vector<IntHalfplane> lines;
    for (auto& hp : h){
        bool same_direction = !lines.empty() && lines.back().cross_dir(hp) == 0;
        if (!same_direction) lines.push_back(hp);
        else if (lines.back().side(hp.px, hp.py) > 0) lines.back() = hp;
    }

    /// The deque is the window dq[head, tail), 13% faster than std::deque at n = 100; the copy only sizes it
    vector<IntHalfplane> dq = lines;
    int head = 0, tail = 0;
    /// A vertex on the new line is popped too, so lines through one point leave one vertex and no zero-length edge
    for (auto& hp : lines){
        while (tail - head > 1 && hp.out(dq[tail - 2], dq[tail - 1])) tail--;
        while (tail - head > 1 && hp.out(dq[head], dq[head + 1])) head++;

        if (tail > head && dq[tail - 1].cross_dir(hp) <= 0) return {};
        dq[tail++] = hp;
    }

    while (tail - head > 2 && dq[head].out(dq[tail - 2], dq[tail - 1])) tail--;
    while (tail - head > 2 && dq[tail - 1].out(dq[head], dq[head + 1])) head++;
    if (tail - head < 3) return {};

    vector<Point> poly;
    for (int i = head; i < tail; i++) poly.push_back(line_intersection(dq[i], dq[i + 1 < tail ? i + 1 : head]));
    return poly;
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
 *   - A float, double or long double read gives strtof's, strtod's or strtold's value, 2-3x faster on short decimals
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

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return 0;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return len;
    }

    /// Digits m and 10^|e| both exact in F round once, as parse does: 2-3x faster than parse on short tokens
    template <typename F, F (*parse)(const char *, char **)>
    bool read_floating(F &x, int max_len, unsigned long long max_mantissa, int max_exp){
        int len = read_token();
        if (!len) return false;

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds float and double twice, long double once
        bool own_width = !FLT_EVAL_METHOD || is_same<F, long double>::value;

        /// A longer token skips the digit scan, which would make its parse about 10% slower
        if (len <= max_len && own_width){
            static const F p10[] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L,
                                    1e10L, 1e11L, 1e12L, 1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L,
                                    1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};
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

            if (!*s && digits && digits <= 19 && m <= max_mantissa && e >= -max_exp && e <= max_exp){
                F v = e < 0 ? m / p10[-e] : m * p10[e];
                x = *tmpbuf == '-' ? -v : v;
                return true;
            }
        }

        x = parse(tmpbuf, nullptr);
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        return read_floating<float, strtof>(x, 10, 1 << 24, 10);
    }

    /// 22 characters keep 16-digit scientific tokens like -1.234567890123456e-05 on the fast path
    bool read_one(double &x){
        return read_floating<double, strtod>(x, 22, 1ULL << 53, 22);
    }

    /// Only an x87 64-bit mantissa takes the fast path: a double-double long double does not round a division once
    bool read_one(long double &x){
        return read_floating<long double, strtold>(x, 21, ~0ULL, LDBL_MANT_DIG == 64 ? 27 : -1);
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

#define MAX 100010
#define BLOCK 2000

int n;
long double X[MAX], bottom[MAX], top[MAX];
vector<Halfplane> active;

vector<Point> solve(int from, int to, vector<Halfplane> &h){
    h = active;
    for (int i = from; i < to; i++){
        h.push_back(Halfplane(Point(0, bottom[i]), Point(1, bottom[i] + X[i])));
        h.push_back(Halfplane(Point(1, top[i] + X[i]), Point(0, top[i])));
    }

    return halfplane_intersection(h, 1e9);
}

bool touches(const Halfplane &hp, const vector<Point> &poly){
    for (const Point &v: poly){
        if (cross(hp.pq, v - hp.p) <= 1e-6 * hp.len) return true;
    }

    return false;
}

int main(){
    int i, done, end, low, high, mid;
    long double max_x = 0, max_y = 0;
    vector<Halfplane> h;

    read(n);
    for (i = 0; i < n; i++){
        read(X[i], bottom[i], top[i]);
        max_x = max(max_x, X[i]), max_y = max(max_y, top[i]);
    }

    for (i = 0; i < n; i++){
        X[i] /= max_x;
        bottom[i] /= (max_y * X[i]), top[i] /= (max_y * X[i]);
    }

    active.push_back(Halfplane(Point(0, 0), Point(0, -1)));
    for (done = 0; done < n; done = end){
        end = min(n, done + BLOCK);
        vector<Point> poly = solve(done, end, h);
        if (poly.empty()) break;

        active.clear();
        for (const Halfplane &hp: h){
            if (touches(hp, poly)) active.push_back(hp);
        }
    }

    if (done < n){
        low = done, high = end - 1;
        while (low < high){
            mid = (low + high + 1) >> 1;
            if (!solve(done, mid, h).empty()) low = mid;
            else high = mid - 1;
        }
        done = low;
    }

    write(done);
    return 0;
}
