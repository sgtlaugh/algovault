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

        int ldo, ldi, rdi, rdo;
        tie(ldo, ldi) = build(l, (l + r) / 2), tie(rdi, rdo) = build((l + r) / 2, r);

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
    int t, n, i, u, v;

    while (read(t)){
        while (t--){
            read(n);
            vector<Point> points(n);
            for (i = 0; i < n; i++) read(points[i].x, points[i].y);

            Delaunay dt(points);
            vector<long long> res(n, LLONG_MAX);
            for (auto& e: dt.edges){
                u = e[0], v = e[1];
                long long dx = points[u].x - points[v].x, dy = points[u].y - points[v].y;
                res[u] = min(res[u], (dx * dx + dy * dy));
                res[v] = min(res[v], (dx * dx + dy * dy));
            }

            for (i = 0; i < n; i++) write(res[i]);
        }
    }

    return 0;
}
