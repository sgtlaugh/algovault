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
        shuffle(order.begin(), order.end(), mt19937(chrono::steady_clock::now().time_since_epoch().count()));

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
 * Fast I/O
 * Buffered stdin and stdout, about 5x faster than scanf for integers, interactive problems included
 *
 * Complexity: O(1) amortized per character, one system call per buffer of up to 8192 bytes
 *
 * fio::flush(): writes the pending output now, which also happens at normal exit
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ..., returns the count read before the input ended
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - Floats print at printf's speed, a floating token read keeps its first 127 characters
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

    bool read_one(double &x){
        if (!read_token()) return false;
        x = strtod(tmpbuf, nullptr);
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

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
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
    int n, i, T = 0;
    char buf[64];

    while (read(n) && n){
        vector<Point3D> points(n);
        for (i = 0; i < n; i++) read(points[i].x, points[i].y, points[i].z);

        ConvexHull3D hull(points);
        snprintf(buf, sizeof(buf), "Case %d: %.2Lf", ++T, hull.area());
        write(string(buf));
    }

    return 0;
}
