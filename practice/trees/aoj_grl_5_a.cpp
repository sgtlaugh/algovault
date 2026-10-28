/***
 *
 * Dynamic Tree Diameter
 * Maintains the diameter of a weighted tree while edge weights change, 0-based nodes
 *
 * Complexity: O(n) to build, O(log n) per weight update, O(1) per diameter query, O(n) memory, about 315 MB at n = 1e6
 *
 * DynamicDiameter tree(n): an empty tree on n >= 1 nodes
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, returns its id, ids are 0, 1, ..., n - 2 in insertion order
 * tree.build(): call after all n - 1 edges
 * tree.update(id, w): sets the weight of edge id to w
 * tree.diameter(): the longest path length, 0 for a single node
 *
 * Notes:
 *   - Weights must be non-negative and their sum must stay at most 2e18
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicDiameter{
    /// g(i): tour depth at i minus the depth before the range; xyz = max g(i) - 2 g(j) + g(k) over i <= j <= k
    struct Node{
        long long sum, x, y, xy, yz, xyz;
    };

    int n, leaves = 0;
    vector<vector<pair<int, int>>> adj;
    vector<long long> weight;
    vector<int> enter, leave;
    vector<Node> seg;

    DynamicDiameter(int n) : n(n), adj(n) {}

    int add_edge(int u, int v, long long w){
        assert(w >= 0);
        int id = weight.size();
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        weight.push_back(w);
        return id;
    }

    void build(){
        assert((int)weight.size() == n - 1);
        enter.assign(n - 1, -1), leave.assign(n - 1, -1);
        vector<int> parent_edge(n, -1), next_edge(n, 0), stack = {0};
        vector<long long> dist(n, 0), tour = {0};

        /// Tour lists a node on entry and again after each child, so edge id covers [enter[id], leave[id] - 1]
        while (!stack.empty()){
            int u = stack.back();
            if (next_edge[u] < (int)adj[u].size()){
                auto [v, id] = adj[u][next_edge[u]++];
                if (id == parent_edge[u]) continue;
                assert(enter[id] == -1);
                parent_edge[v] = id, dist[v] = dist[u] + weight[id];
                enter[id] = tour.size(), tour.push_back(dist[v]);
                stack.push_back(v);
                continue;
            }
            stack.pop_back();
            if (stack.empty()) break;
            leave[parent_edge[u]] = tour.size(), tour.push_back(dist[stack.back()]);
        }
        assert((int)tour.size() == 2 * n - 1);

        /// A power of two keeps every leaf on one level, so seg[1] combines the whole tour left to right
        leaves = 1;
        while (leaves < (int)tour.size()) leaves *= 2;
        seg.assign(2 * leaves, Node{});
        for (int i = 1; i < (int)tour.size(); i++) seg[leaves + i] = leaf(tour[i] - tour[i - 1]);
        for (int i = leaves - 1; i >= 1; i--) seg[i] = combine(seg[2 * i], seg[2 * i + 1]);
    }

    /// Both leaves change before any parent is recombined, so no node ever sees a half-applied update
    void update(int id, long long w){
        assert(w >= 0);
        long long delta = w - weight[id];
        weight[id] = w;
        int l = leaves + enter[id], r = leaves + leave[id];
        seg[l] = leaf(seg[l].sum + delta), seg[r] = leaf(seg[r].sum - delta);
        while (l > 1){
            l /= 2, r /= 2;
            seg[l] = combine(seg[2 * l], seg[2 * l + 1]);
            if (r != l) seg[r] = combine(seg[2 * r], seg[2 * r + 1]);
        }
    }

    long long diameter() const{
        return seg[1].xyz;
    }

    static Node leaf(long long d){
        return {d, d, -2 * d, -d, -d, 0};
    }

    /// r's values are relative to the end of l, so each term shifts by l.sum once per g it holds
    static Node combine(const Node& l, const Node& r){
        long long s = l.sum;
        Node res;
        res.sum = s + r.sum;
        res.x = max(l.x, r.x + s);
        res.y = max(l.y, r.y - 2 * s);
        res.xy = max({l.xy, r.xy - s, l.x + (r.y - 2 * s)});
        res.yz = max({l.yz, r.yz - s, l.y + (r.x + s)});
        res.xyz = max({l.xyz, r.xyz, l.xy + (r.x + s), l.x + (r.yz - s)});
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
    int n, i, u, v, w;

    read(n);
    DynamicDiameter tree(n);
    for (i = 1; i < n; i++){
        read(u, v, w);
        tree.add_edge(u, v, w);
    }

    tree.build();
    write(tree.diameter());
    return 0;
}
