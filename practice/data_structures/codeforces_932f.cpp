/***
 *
 * Li Chao Tree
 * Holds lines y = k * x + b over integer x in [lo, hi] and returns the minimum or maximum at a point
 *
 * Complexity:
 *   - LiChaoTree: O(log C) per add_line and query, O(log^2 C) per add_segment, C = hi - lo + 1
 *   - PersistentLiChaoTree: O(log C) per add_line and query, add_line makes at most log2(C) + 2 nodes
 *
 * LiChaoTree<T> tree(lo, hi): no lines yet, x in [lo, hi], the range can span up to [-1e18, 1e18]
 * LiChaoTree<T, true> tree(lo, hi): the same tree answering the maximum
 * tree.add_line(k, b): add y = k * x + b for every x in [lo, hi]
 * tree.add_segment(k, b, l, r): add the line only for x in [l, r] inclusive, the part outside [lo, hi] is ignored
 * tree.query(x): minimum over the lines covering x, or NONE if no line covers it, lo <= x <= hi
 *
 * PersistentLiChaoTree<T> pt(lo, hi): version 0 holds no lines, same range limits as LiChaoTree
 * PersistentLiChaoTree<T, true> pt(lo, hi): the same tree answering the maximum
 * pt.add_line(version, k, b): a new version with the lines of version plus y = k * x + b, version stays intact
 * pt.query(version, x): minimum over the lines of version at x, or NONE if it has no line, lo <= x <= hi
 *
 * Notes:
 *   - NONE is the largest value of T, or the lowest when answering the maximum
 *   - Every k * x + b over the range must fit in T, e.g. |k|, |x| <= 1e9 and |b| <= 1e18 for long long
 *   - LiChaoTree makes nodes on demand, at most 1 per add_line and O(log C) per add_segment
 *   - A node is 24 bytes for long long, random segments at C = 2e9 average ~31 nodes, so 4e5 of them take ~400 MB peak
 *   - Convex hull trick: for dp[i] = min over j < i of dp[j] + a[j] * x[i], query(x[i]) then add_line(a[i], dp[i])
 *   - On root to node paths: dp[v] = pt.query(ver[parent], x[v]), then ver[v] = pt.add_line(ver[parent], a[v], dp[v])
 *   - The persistent DFS needs no rollback when it backtracks, siblings branch off the same parent version
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, bool MAXIMIZE = false>
struct LiChaoTree{
    static constexpr T NONE = MAXIMIZE ? numeric_limits<T>::lowest() : numeric_limits<T>::max();

    struct Line{
        T k, b;

        T eval(long long x) const{
            return k * x + b;
        }
    };

    /// A node on a segment's split path holds y = NONE, which never wins, so no node needs a has-line flag
    struct Node{
        Line line;
        int left, right;
    };

    long long lo, hi;
    int root = -1;
    vector<Node> nodes;

    LiChaoTree(long long lo, long long hi) : lo(lo), hi(hi) {}

    void add_line(T k, T b){
        root = insert(root, lo, hi, {k, b});
    }

    void add_segment(T k, T b, long long l, long long r){
        if (l <= r) root = insert_segment(root, lo, hi, l, r, {k, b});
    }

    T query(long long x) const{
        assert(lo <= x && x <= hi);

        T res = NONE;
        long long l = lo, r = hi;
        for (int node = root; node != -1;){
            const Node& cur = nodes[node];
            T y = cur.line.eval(x);
            if (better(y, res)) res = y;

            long long m = l + (r - l) / 2;
            if (x <= m) node = cur.left, r = m;
            else node = cur.right, l = m + 1;
        }

        return res;
    }

    /// The early stop needs functions that cross at most once
    int insert(int node, long long l, long long r, Line line){
        if (node == -1) return new_node(line);

        for (int top = node; ; ){
            Line& cur = nodes[node].line;
            bool left_better = better(line.eval(l), cur.eval(l));
            bool right_better = better(line.eval(r), cur.eval(r));
            if (left_better == right_better){
                if (left_better) cur = line;
                return top;
            }

            long long m = l + (r - l) / 2;
            bool mid_better = better(line.eval(m), cur.eval(m));
            if (mid_better) swap(cur, line);

            bool go_left = left_better != mid_better;
            int child = go_left ? nodes[node].left : nodes[node].right;
            if (child == -1){
                child = new_node(line);
                (go_left ? nodes[node].left : nodes[node].right) = child;
                return top;
            }

            node = child;
            if (go_left) r = m;
            else l = m + 1;
        }
    }

    int insert_segment(int node, long long l, long long r, long long ql, long long qr, const Line& line){
        if (qr < l || r < ql) return node;
        if (ql <= l && r <= qr) return insert(node, l, r, line);
        if (node == -1) node = new_node({0, NONE});

        long long m = l + (r - l) / 2;
        int left = insert_segment(nodes[node].left, l, m, ql, qr, line);  /// can reallocate nodes, assign after
        nodes[node].left = left;
        int right = insert_segment(nodes[node].right, m + 1, r, ql, qr, line);
        nodes[node].right = right;
        return node;
    }

    static bool better(T a, T b){
        return MAXIMIZE ? a > b : a < b;
    }

    int new_node(const Line& line){
        nodes.push_back({line, -1, -1});
        return (int)nodes.size() - 1;
    }
};

template <typename T, bool MAXIMIZE = false>
struct PersistentLiChaoTree{
    using Base = LiChaoTree<T, MAXIMIZE>;
    using Line = typename Base::Line;
    using Node = typename Base::Node;

    static constexpr T NONE = Base::NONE;

    long long lo, hi;
    vector<int> roots = {-1};
    vector<Node> nodes;

    PersistentLiChaoTree(long long lo, long long hi) : lo(lo), hi(hi) {}

    /// Midpoint rule of LiChaoTree::insert without its early stop, writing into a copy of every node it visits
    int add_line(int version, T k, T b){
        assert(0 <= version && version < (int)roots.size());

        Line line = {k, b};
        int node = roots[version], top = -1, parent = -1;
        bool from_left = false;
        long long l = lo, r = hi;

        while (true){
            Node copy = node == -1 ? Node{line, -1, -1} : nodes[node];
            nodes.push_back(copy);
            int cur = (int)nodes.size() - 1;
            if (parent == -1) top = cur;
            else (from_left ? nodes[parent].left : nodes[parent].right) = cur;
            if (node == -1) break;

            long long m = l + (r - l) / 2;
            Line& held = nodes[cur].line;
            bool left_better = Base::better(line.eval(l), held.eval(l));
            bool mid_better = Base::better(line.eval(m), held.eval(m));
            if (mid_better) swap(held, line);
            if (l == r) break;

            parent = cur, from_left = left_better != mid_better;
            node = from_left ? nodes[cur].left : nodes[cur].right;
            if (from_left) r = m;
            else l = m + 1;
        }

        roots.push_back(top);
        return (int)roots.size() - 1;
    }

    T query(int version, long long x) const{
        assert(0 <= version && version < (int)roots.size());
        assert(lo <= x && x <= hi);

        T res = NONE;
        long long l = lo, r = hi;
        for (int node = roots[version]; node != -1;){
            const Node& cur = nodes[node];
            T y = cur.line.eval(x);
            if (Base::better(y, res)) res = y;

            long long m = l + (r - l) / 2;
            if (x <= m) node = cur.left, r = m;
            else node = cur.right, l = m + 1;
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

#define MAX 100010

int n, a[MAX], b[MAX], parent[MAX], id[MAX];
long long dp[MAX];
vector<int> adj[MAX];
vector<pair<long long, long long>> lines[MAX];
vector<LiChaoTree<long long>> trees;

vector<int> bfs(){
    vector<int> order = {1};

    for (int i = 0; i < (int)order.size(); i++){
        int u = order[i];
        for (int v: adj[u]){
            if (v != parent[u]) parent[v] = u, order.push_back(v);
        }
    }

    return order;
}

void merge(int u, int v){
    if (lines[id[u]].size() < lines[id[v]].size()) swap(id[u], id[v]);

    int x = id[u], y = id[v];
    for (auto line: lines[y]){
        trees[x].add_line(line.first, line.second);
        lines[x].push_back(line);
    }

    vector<pair<long long, long long>>().swap(lines[y]);
    vector<LiChaoTree<long long>::Node>().swap(trees[y].nodes);
}

int main(){
    int i, u, v;

    read(n);
    for (i = 1; i <= n; i++) read(a[i]);
    for (i = 1; i <= n; i++) read(b[i]);
    for (i = 1; i < n; i++){
        read(u, v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> order = bfs();
    trees.assign(n + 1, LiChaoTree<long long>(-100000, 100000));
    for (i = 1; i <= n; i++) id[i] = i;

    for (i = n - 1; i >= 0; i--){
        u = order[i];
        if (lines[id[u]].empty()) dp[u] = 0;
        else dp[u] = trees[id[u]].query(a[u]);

        trees[id[u]].add_line(b[u], dp[u]);
        lines[id[u]].push_back({b[u], dp[u]});
        if (u != 1) merge(parent[u], u);
    }

    vector<long long> res(dp + 1, dp + n + 1);
    write(res);
    return 0;
}
