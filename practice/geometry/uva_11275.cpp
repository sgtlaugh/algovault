/***
 *
 * 3D Geometry
 * Distances between points, lines, segments, planes and triangles in space, line-plane intersection and volume
 *
 * Complexity: O(1) per call, O(number of face vertices) for polyhedron_volume
 *
 * Point3 of doubles with +, -, * and / by a scalar; a plane or triangle is given by three points a, b, c
 * cross(u, v): the cross product u x v
 * dist(p, q): Euclidean distance between two points
 * dist_point_line(p, a, b): distance from p to the line through a != b
 * dist_point_plane(p, a, b, c): signed distance, positive on the side of cross(b - a, c - a), a, b, c not collinear
 * dist_point_segment(p, a, b): distance from p to segment ab, a == b allowed
 * dist_point_triangle(p, a, b, c): distance from p to the filled triangle abc, degenerate triangles allowed
 * dist_segment_segment(a, b, c, d): distance between segments ab and cd, zero-length ones allowed
 * dist_segment_triangle(p, q, a, b, c): distance between segment pq and the filled triangle abc, about 0 when they meet
 * dot(u, v): the dot product
 * line_plane_intersection(p, q, a, b, c): {1, point} for one point, {0, p} parallel, {-1, p} the line lies in the plane
 * norm(u): length of u
 * polyhedron_volume(points, faces): volume of a closed polyhedron, each face counter-clockwise seen from outside
 *
 * Notes:
 *   - Every result is a double; distances hold to 1e-7 at coordinates up to 1e6, nearly parallel segments included
 *   - A line is parallel to the plane when the sine of their angle is below EPS = 1e-9
 *   - Two triangles share a point iff some edge of one is within the tolerance of the other: dist_segment_triangle
 *   - For a convex hull of integer points use convex_hull_3D.cpp, whose volume() is exact
 *   - Defines EPS = 1e-9 as geometry.cpp does: pasted next to it, keep one EPS
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const double EPS = 1e-9;

struct Point3{
    double x, y, z;

    Point3 operator+(const Point3& p) const{
        return {x + p.x, y + p.y, z + p.z};
    }

    Point3 operator-(const Point3& p) const{
        return {x - p.x, y - p.y, z - p.z};
    }

    Point3 operator*(double k) const{
        return {x * k, y * k, z * k};
    }

    Point3 operator/(double k) const{
        return {x / k, y / k, z / k};
    }
};

Point3 cross(const Point3& u, const Point3& v){
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

double dot(const Point3& u, const Point3& v){
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

double norm(const Point3& u){
    return sqrt(dot(u, u));
}

double dist(const Point3& p, const Point3& q){
    return norm(p - q);
}

double dist_point_line(const Point3& p, const Point3& a, const Point3& b){
    return norm(cross(b - a, p - a)) / norm(b - a);
}

double dist_point_plane(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    return dot(n, p - a) / norm(n);
}

double dist_point_segment(const Point3& p, const Point3& a, const Point3& b){
    Point3 ab = b - a;
    double len2 = dot(ab, ab), t = len2 > 0 ? dot(p - a, ab) / len2 : 0;
    return dist(p, a + ab * max(0.0, min(1.0, t)));
}

/// No tolerance needed: a point misjudged near an edge gets that edge's distance from the caller anyway
bool projects_inside(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    if (dot(n, n) == 0) return false;
    return dot(cross(b - a, p - a), n) >= 0 && dot(cross(c - b, p - b), n) >= 0 && dot(cross(a - c, p - c), n) >= 0;
}

double dist_point_triangle(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    if (projects_inside(p, a, b, c)) return fabs(dist_point_plane(p, a, b, c));
    return min({dist_point_segment(p, a, b), dist_point_segment(p, b, c), dist_point_segment(p, c, a)});
}

double dist_segment_segment(const Point3& a, const Point3& b, const Point3& c, const Point3& d){
    Point3 u = b - a, v = d - c, r = a - c;
    double uu = dot(u, u), vv = dot(v, v), uv = dot(u, v), ur = dot(u, r), vr = dot(v, r);
    if (uu == 0) return dist_point_segment(a, c, d);
    if (vv == 0) return dist_point_segment(c, a, b);

    /// (u x v).(v x r) equals uv vr - vv ur without its cancellation, which cost 1e-6 on long near-parallel segments
    Point3 w = cross(u, v);
    double den = dot(w, w), s = den > 0 ? max(0.0, min(1.0, dot(w, cross(v, r)) / den)) : 0;
    double t = (uv * s + vr) / vv;
    if (t < 0) t = 0, s = max(0.0, min(1.0, -ur / uu));
    else if (t > 1) t = 1, s = max(0.0, min(1.0, (uv - ur) / uu));

    return dist(a + u * s, c + v * t);
}

double dist_segment_triangle(const Point3& p, const Point3& q, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    double sp = dot(n, p - a), sq = dot(n, q - a);

    if ((sp < 0 && sq > 0) || (sp > 0 && sq < 0)){
        if (projects_inside(p + (q - p) * (sp / (sp - sq)), a, b, c)) return 0;
    }

    double res = min(dist_point_triangle(p, a, b, c), dist_point_triangle(q, a, b, c));
    return min({res, dist_segment_segment(p, q, a, b), dist_segment_segment(p, q, b, c), dist_segment_segment(p, q, c, a)});
}

pair<int, Point3> line_plane_intersection(const Point3& p, const Point3& q, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a), d = q - p;
    double den = dot(n, d);
    if (fabs(den) > EPS * norm(n) * norm(d)) return {1, p + d * (dot(n, a - p) / den)};

    bool in_plane = fabs(dist_point_plane(p, a, b, c)) <= EPS * max(norm(p - a), norm(d));
    return {in_plane ? -1 : 0, p};
}

double polyhedron_volume(const vector<Point3>& points, const vector<vector<int>>& faces){
    double res = 0;
    for (const vector<int>& f : faces){
        /// Fanned from the face's first vertex: each term is 6 times a signed tetrahedron against the origin
        for (int i = 1; i + 1 < (int)f.size(); i++) res += dot(points[f[0]], cross(points[f[i]], points[f[i + 1]]));
    }

    return res / 6;
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
        while (isdigit(c)){
            v = v * 10 - (c - '0');

            /// The digits already buffered are read through a local pointer, 15% faster than a read_char call per digit
            const char *p = inbuf + inptr, *end = inbuf + buf_len;
            while (p < end && (unsigned)(*p - '0') < 10) v = v * 10 - (*p++ - '0');
            inptr = p - inbuf;

            /// A buffered non-digit ends the number and is eaten, else the digits may go on after a refill
            if (p < end){
                inptr++;
                break;
            }
            c = read_char();
        }
        x = neg ? v : -v;

        return true;
    }

    /// The C locale's isspace set, false for EOF, and 1.4-1.7x faster for strings and chars than isspace's table lookup
    inline bool is_space(int c){
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (is_space(c)) c = read_char();
        if (c == EOF) return 0;

        for (; !is_space(c) && c != EOF; c = read_char()){
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

        while (is_space(c)) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !is_space(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (is_space(c));
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

Point3 tri[2][3];

bool touch(){
    int i, j;

    for (i = 0; i < 2; i++){
        for (j = 0; j < 3; j++){
            const Point3* t = tri[i ^ 1];
            if (dist_segment_triangle(tri[i][j], tri[i][(j + 1) % 3], t[0], t[1], t[2]) <= 1e-6 + 1e-12) return true;
        }
    }

    return false;
}

int main(){
    int t, i, j;

    while (read(t)){
        while (t--){
            for (i = 0; i < 2; i++){
                for (j = 0; j < 3; j++) read(tri[i][j].x, tri[i][j].y, tri[i][j].z);
            }

            write(touch() ? 1 : 0);
        }
    }

    return 0;
}
