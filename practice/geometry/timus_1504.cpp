/***
 *
 * Delaunay Triangulation
 * Triangulates points so that no point lies strictly inside any triangle's circumcircle
 *
 * Complexity: O(n log n)
 *
 * Delaunay dt(points): the Delaunay triangulation of points in any order, duplicates allowed
 * dt.triangles: index triples {a, b, c} into points, counter-clockwise with a the smallest
 * dt.edges: index pairs {u, v}, every edge of the triangulation once
 *
 * Notes:
 *   - |x|, |y| <= 1e9, requires __int128
 *   - Collinear points give no triangles, and edges join consecutive points along the line
 *   - Four or more cocircular points get any one of their valid triangulations
 *   - Duplicates: triangles and edges use the copy with the smallest index, edges also join every other copy u to it
 *   - edges contains a Euclidean MST: Kruskal from minimum_spanning_tree.cpp with weight dist2(u, v) finds it
 *   - Defines its own Point: pasted next to another geometry file, merge the Point blocks first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;
};

/// Guibas-Stolfi divide and conquer on quad-edges, adapted from KACTL FastDelaunay.h (MIT licence)
struct Delaunay{
    vector<Point> pts;
    vector<int> id;
    /// Directed edge e runs org[e] -> org[e ^ 1], onext and oprev link the edges leaving org[e] in ccw order
    vector<int> org, onext, oprev, free_edges;
    vector<array<int, 3>> triangles;
    vector<array<int, 2>> edges;

    Delaunay(const vector<Point>& points){
        int n = points.size();
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){
            return make_tuple(points[a].x, points[a].y, a) < make_tuple(points[b].x, points[b].y, b);
        });

        for (int i : order){
            const Point& p = points[i];
            if (!pts.empty() && pts.back().x == p.x && pts.back().y == p.y) edges.push_back({i, id.back()});
            else pts.push_back(p), id.push_back(i);
        }

        if (pts.size() >= 2) build(0, pts.size());
        collect();
    }

    long long cross(int o, int a, int b) const{
        return (pts[a].x - pts[o].x) * (pts[b].y - pts[o].y) - (pts[a].y - pts[o].y) * (pts[b].x - pts[o].x);
    }

    /// True if d is strictly inside the circle through ccw a, b, c; partial sums reach 1.3e38 < 2^127 at |x|, |y| = 1e9
    bool in_circle(int a, int b, int c, int d) const{
        auto lift = [&](int u){
            long long dx = pts[u].x - pts[d].x, dy = pts[u].y - pts[d].y;
            return (__int128)(dx * dx + dy * dy);
        };
        return lift(a) * cross(d, b, c) + lift(b) * cross(d, c, a) + lift(c) * cross(d, a, b) > 0;
    }

    int dest(int e) const{
        return org[e ^ 1];
    }

    int lnext(int e) const{
        return oprev[e ^ 1];
    }

    /// Returns {ccw hull edge out of the leftmost point, cw hull edge out of the rightmost point} of pts[l, r)
    pair<int, int> build(int l, int r){
        if (r - l == 2){
            int a = make_edge(l, l + 1);
            return {a, a ^ 1};
        }
        if (r - l == 3){
            int a = make_edge(l, l + 1), b = make_edge(l + 1, l + 2);
            splice(a ^ 1, b);
            long long turn = cross(l, l + 1, l + 2);
            if (turn == 0) return {a, b ^ 1};
            int c = connect(b, a);
            if (turn > 0) return {a, b ^ 1};
            return {c ^ 1, c};
        }

        auto [ldo, ldi] = build(l, (l + r) / 2);
        auto [rdi, rdo] = build((l + r) / 2, r);

        while (true){
            if (cross(org[ldi], dest(ldi), org[rdi]) > 0) ldi = lnext(ldi);
            else if (cross(org[rdi], dest(rdi), org[ldi]) < 0) rdi = onext[rdi ^ 1];
            else break;
        }

        int base = connect(rdi ^ 1, ldi);
        if (org[ldi] == org[ldo]) ldo = base ^ 1;
        if (org[rdi] == org[rdo]) rdo = base;

        /// Zips the two halves bottom up, base runs from the right half to the left half along the lower boundary
        auto valid = [&](int e){ return cross(dest(e), dest(base), org[base]) > 0; };
        while (true){
            int lcand = onext[base ^ 1];
            while (valid(lcand) && in_circle(dest(base), org[base], dest(lcand), dest(onext[lcand]))){
                int next = onext[lcand];
                remove_edge(lcand);
                lcand = next;
            }

            int rcand = oprev[base];
            while (valid(rcand) && in_circle(dest(base), org[base], dest(rcand), dest(oprev[rcand]))){
                int next = oprev[rcand];
                remove_edge(rcand);
                rcand = next;
            }

            if (!valid(lcand) && !valid(rcand)) break;
            if (!valid(lcand) || (valid(rcand) && in_circle(dest(lcand), org[lcand], org[rcand], dest(rcand)))){
                base = connect(rcand, base ^ 1);
            }
            else base = connect(base ^ 1, lcand ^ 1);
        }

        return {ldo, rdo};
    }

    int make_edge(int u, int v){
        int e = org.size();
        if (!free_edges.empty()){
            e = free_edges.back();
            free_edges.pop_back();
        }
        else{
            org.resize(e + 2);
            onext.resize(e + 2);
            oprev.resize(e + 2);
        }

        org[e] = u, org[e ^ 1] = v;
        onext[e] = oprev[e] = e;
        onext[e ^ 1] = oprev[e ^ 1] = e ^ 1;
        return e;
    }

    /// Merges the rings around org[a] and org[b] if they differ, else splits them, as in Guibas-Stolfi
    void splice(int a, int b){
        int an = onext[a], bn = onext[b];
        onext[a] = bn, onext[b] = an;
        oprev[bn] = a, oprev[an] = b;
    }

    /// Adds the edge dest(a) -> org[b], with a's left face to its left and b's left face to its right
    int connect(int a, int b){
        int e = make_edge(dest(a), org[b]);
        splice(e, lnext(a));
        splice(e ^ 1, b);
        return e;
    }

    void remove_edge(int e){
        splice(e, oprev[e]);
        splice(e ^ 1, oprev[e ^ 1]);
        free_edges.push_back(e & ~1);
    }

    void collect(){
        int m = org.size();
        vector<bool> alive(m, true);
        for (int e : free_edges) alive[e] = alive[e ^ 1] = false;

        for (int e = 0; e < m; e += 2){
            if (alive[e]) edges.push_back({id[org[e]], id[org[e ^ 1]]});
        }

        /// Each triangle is a ccw face of three edges, reported once from its smallest edge; the outer face turns cw
        for (int a = 0; a < m; a++){
            if (!alive[a]) continue;
            int b = lnext(a), c = lnext(b);
            if (lnext(c) == a && a < b && a < c && cross(org[a], org[b], org[c]) > 0){
                array<int, 3> t = {id[org[a]], id[org[b]], id[org[c]]};
                rotate(t.begin(), min_element(t.begin(), t.end()), t.end());
                triangles.push_back(t);
            }
        }
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

const long double SCALE = 1e6, EPS = 1e-3;

int n, w[1010];
long double R;
vector<Point> pts;
vector<long double> angles;

int value(long double a){
    long double x = R * cosl(a), y = R * sinl(a), best = 1e30;
    int i, res = 0;

    vector<long double> d(n);
    for (i = 0; i < n; i++){
        d[i] = hypotl(x - pts[i].x, y - pts[i].y);
        best = min(best, d[i]);
    }
    for (i = 0; i < n; i++){
        if (d[i] <= best + EPS) res = max(res, w[i]);
    }

    return res;
}

void bisector(Point u, Point v){
    long double nx = v.x - u.x, ny = v.y - u.y, len = hypotl(nx, ny);
    long double c = ((long double)v.x * v.x + (long double)v.y * v.y - (long double)u.x * u.x - (long double)u.y * u.y) / 2 / len;

    if (fabsl(c) > R) return;
    long double h = sqrtl(R * R - c * c);
    nx /= len, ny /= len;
    angles.push_back(atan2l(ny * c + nx * h, nx * c - ny * h));
    angles.push_back(atan2l(ny * c - nx * h, nx * c + ny * h));
}

int main(){
    int i, k, cur, best;
    long double x, y, a, res;

    read(R, n);
    R *= SCALE;
    pts.resize(n);
    for (i = 0; i < n; i++){
        read(x, y, w[i]);
        pts[i] = {llroundl(x * SCALE), llroundl(y * SCALE)};
    }

    Delaunay dt(pts);
    for (auto e: dt.edges) bisector(pts[e[0]], pts[e[1]]);
    sort(angles.begin(), angles.end());

    k = angles.size();
    res = 0, best = value(0);
    for (i = 0; i < k; i++){
        a = (i + 1) < k ? (angles[i] + angles[i + 1]) / 2 : (angles[i] + angles[0]) / 2 + M_PI;
        cur = value(a);
        if (cur > best) best = cur, res = a;
    }
    for (i = 0; i < k; i++){
        cur = value(angles[i]);
        if (cur > best) best = cur, res = angles[i];
    }

    write(fixed(R * cosl(res) / SCALE, 9), fixed(R * sinl(res) / SCALE, 9));
    return 0;
}
