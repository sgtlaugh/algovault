/***
 *
 * KD-Tree
 * Static 2D point set: the point nearest to a query by Euclidean or Manhattan distance, and the points in a rectangle
 *
 * Complexity:
 *   - KDTree: O(n log n)
 *   - nearest: about O(log n) on random or clustered points, O(n) for a query at the centre of points on a circle
 *   - report: about O(sqrt n + k) for k reported points
 *
 * KDTree tree(points): tree over points whose nearest uses the squared Euclidean distance dx^2 + dy^2
 * KDTree tree(points, KDTree::MANHATTAN): tree over points whose nearest uses the Manhattan distance |dx| + |dy|
 * tree.nearest(q): {dist, i}, a point points[i] nearest to q and its distance, {LLONG_MAX, -1} if there is none
 * tree.nearest(q, skip): {dist, i}, a point nearest to q among all points except points[skip], i = -1 if there is none
 * tree.report(x1, y1, x2, y2): indices of the points with x1 <= x <= x2 and y1 <= y <= y2, in no particular order
 *
 * Notes:
 *   - |x|, |y| <= 1e9 for points and queries, so a squared distance (at most 8e18) fits in a long long
 *   - nearest(points[i], i) is the nearest other point, at dist 0 when points[i] has a duplicate
 *   - Ties between equally near points are broken arbitrarily
 *   - Reimplemented from the idea of the KACTL KD-tree (no licence stated): bounding box pruning, nearer child first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Point{
    long long x, y;
};

struct KDTree{
    enum Metric{EUCLIDEAN, MANHATTAN};

    struct Neighbor{
        long long dist;
        int i;
    };

    /// The subtree over positions [l, r) is rooted at (l + r) / 2, which holds the box of the whole subtree
    struct Node{
        Point p;
        int i;
        long long x1, x2, y1, y2;
    };

    Metric metric;
    vector<Node> nodes;

    KDTree(const vector<Point>& points, Metric metric = EUCLIDEAN) : metric(metric), nodes(points.size()){
        for (int i = 0; i < (int)points.size(); i++) nodes[i].p = points[i], nodes[i].i = i;
        build(0, nodes.size());
    }

    Neighbor nearest(const Point& q, int skip = -1) const {
        Neighbor best = {LLONG_MAX, -1};
        if (!nodes.empty()) search(0, nodes.size(), q, skip, best);
        return best;
    }

    vector<int> report(long long x1, long long y1, long long x2, long long y2) const {
        vector<int> res;
        collect(0, nodes.size(), x1, y1, x2, y2, res);
        return res;
    }

    void build(int l, int r){
        if (l >= r) return;

        long long x1 = LLONG_MAX, x2 = LLONG_MIN, y1 = LLONG_MAX, y2 = LLONG_MIN;
        for (int k = l; k < r; k++){
            x1 = min(x1, nodes[k].p.x), x2 = max(x2, nodes[k].p.x);
            y1 = min(y1, nodes[k].p.y), y2 = max(y2, nodes[k].p.y);
        }

        int m = (l + r) / 2;
        bool by_x = x2 - x1 >= y2 - y1;
        nth_element(nodes.begin() + l, nodes.begin() + m, nodes.begin() + r, [&](const Node& a, const Node& b){
            return by_x ? a.p.x < b.p.x : a.p.y < b.p.y;
        });
        nodes[m].x1 = x1, nodes[m].x2 = x2, nodes[m].y1 = y1, nodes[m].y2 = y2;

        build(l, m);
        build(m + 1, r);
    }

    void search(int l, int r, const Point& q, int skip, Neighbor& best) const {
        int m = (l + r) / 2;
        const Node& node = nodes[m];
        if (node.i != skip){
            long long d = dist(node.p.x - q.x, node.p.y - q.y);
            if (d < best.dist) best = {d, node.i};
        }

        long long dl = box_dist(l, m, q), dr = box_dist(m + 1, r, q);
        if (dl <= dr){
            if (dl < best.dist) search(l, m, q, skip, best);
            if (dr < best.dist) search(m + 1, r, q, skip, best);
        }
        else{
            if (dr < best.dist) search(m + 1, r, q, skip, best);
            if (dl < best.dist) search(l, m, q, skip, best);
        }
    }

    /// LLONG_MAX for an empty range, which no best distance can beat, so empty children are never visited
    long long box_dist(int l, int r, const Point& q) const {
        if (l >= r) return LLONG_MAX;

        const Node& node = nodes[(l + r) / 2];
        long long dx = max({node.x1 - q.x, q.x - node.x2, 0LL});
        long long dy = max({node.y1 - q.y, q.y - node.y2, 0LL});
        return dist(dx, dy);
    }

    long long dist(long long dx, long long dy) const {
        return metric == EUCLIDEAN ? dx * dx + dy * dy : abs(dx) + abs(dy);
    }

    void collect(int l, int r, long long x1, long long y1, long long x2, long long y2, vector<int>& res) const {
        if (l >= r) return;

        int m = (l + r) / 2;
        const Node& node = nodes[m];
        if (node.x2 < x1 || x2 < node.x1 || node.y2 < y1 || y2 < node.y1) return;
        if (x1 <= node.x1 && node.x2 <= x2 && y1 <= node.y1 && node.y2 <= y2){
            for (int k = l; k < r; k++) res.push_back(nodes[k].i);
            return;
        }

        if (x1 <= node.p.x && node.p.x <= x2 && y1 <= node.p.y && node.p.y <= y2) res.push_back(node.i);
        collect(l, m, x1, y1, x2, y2, res);
        collect(m + 1, r, x1, y1, x2, y2, res);
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
    int t, n, i;

    while (read(t)){
        while (t--){
            read(n);
            vector<Point> points(n);
            for (i = 0; i < n; i++) read(points[i].x, points[i].y);

            KDTree tree(points);
            vector<long long> res(n);
            for (auto& node: tree.nodes) res[node.i] = tree.nearest(node.p, node.i).dist;
            for (i = 0; i < n; i++) write(res[i]);
        }
    }

    return 0;
}
