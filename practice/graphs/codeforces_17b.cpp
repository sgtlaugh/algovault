/***
 *
 * Directed MST
 * Constructs a rooted tree of minimum total weight from the root node, with the parent of every node
 *
 * Complexity: O(m log m)
 *
 * DirectedMST g(n): nodes are 0 to n - 1
 * g.add_edge(u, v, w): directed edge u -> v with long long weight w, negative allowed
 * g.solve(root): {exists, weight, parent}, parent[v] the tail of v's tree edge and -1 at root
 *
 * Notes:
 *   - Requires n max|w| <= 1e18
 *   - exists is false when some node is unreachable from root, then weight is 0 and parent empty
 *   - solve is const, so one graph answers several roots
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Arborescence{
    bool exists;
    long long weight;
    vector<int> parent;
};

struct DirectedMST{
    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;

    DirectedMST(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    Arborescence solve(int root) const{
        struct Node{
            long long key, lazy;
            int l, r;
        };

        int m = edges.size();
        vector<Node> t(m);
        vector<int> heap(2 * n, -1), dsu(2 * n), up(2 * n, -1), in(2 * n, -1), seen(2 * n, -1), path(n), spine;
        iota(dsu.begin(), dsu.end(), 0);

        auto find = [&](int x){
            while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]];
            return x;
        };

        auto push = [&](int a){
            if (!t[a].lazy) return;  /// most nodes carry no pending add, skipping saves two scattered child writes
            if (t[a].l != -1) t[t[a].l].lazy += t[a].lazy;
            if (t[a].r != -1) t[t[a].r].lazy += t[a].lazy;
            t[a].key += t[a].lazy, t[a].lazy = 0;
        };

        /// Iterative because a single skew heap merge can walk an O(m) right spine
        auto merge = [&](int a, int b){
            spine.clear();
            while (a != -1 && b != -1){
                push(a), push(b);
                if (t[a].key > t[b].key) swap(a, b);
                spine.push_back(a);
                a = t[a].r;
            }

            int res = a != -1 ? a : b;
            while (!spine.empty()){
                int x = spine.back();
                spine.pop_back();
                t[x].r = t[x].l, t[x].l = res, res = x;
            }
            return res;
        };

        for (int i = 0; i < m; i++){
            t[i] = {edges[i].w, 0, -1, -1};
            heap[edges[i].v] = merge(heap[edges[i].v], i);
        }

        long long weight = 0;
        int nodes = n;
        seen[root] = root;
        for (int s = 0; s < n; s++){
            int u = s, len = 0;
            while (seen[u] == -1){
                if (heap[u] == -1) return {false, 0, {}};

                int e = heap[u];
                push(e);
                weight += t[e].key;
                heap[u] = merge(t[e].l, t[e].r);
                if (heap[u] != -1) t[heap[u]].lazy -= t[e].key;

                in[u] = e, path[len++] = u, seen[u] = s;
                u = find(edges[e].u);
                if (seen[u] != s) continue;

                /// A self loop of a contracted node is no cycle, popping it again just takes its next edge
                if (u == path[len - 1]){
                    len--, seen[u] = -1;
                    continue;
                }

                int cycle = nodes++, x;
                do{
                    x = path[--len];
                    heap[cycle] = merge(heap[cycle], heap[x]);
                    dsu[x] = up[x] = cycle;
                } while (x != u);
                u = cycle;
            }
        }

        /// Parents are built before children, so a node is known to be dropped before it is reached
        vector<int> parent(n, -1);
        vector<bool> dropped(nodes, false);
        for (int x = nodes - 1; x >= 0; x--){
            if (x == root || dropped[x]) continue;
            int e = in[x], v = edges[e].v;
            parent[v] = edges[e].u;
            for (int y = v; y != x; y = up[y]) dropped[y] = true;
        }

        return {true, weight, parent};
    }
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
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
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
            if (buf_len == 0) return EOF;
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

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

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

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
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

#define MAX 1010

int n, q[MAX];

int main(){
    int m, i, a, b, c, root = 0, cnt = 0;

    read(n);
    for (i = 0; i < n; i++){
        read(q[i]);
        if (q[i] > q[root]) root = i;
    }
    for (i = 0; i < n; i++){
        if (q[i] == q[root]) cnt++;
    }

    read(m);
    DirectedMST g(n);
    for (i = 0; i < m; i++){
        read(a, b, c);
        g.add_edge(a - 1, b - 1, c);
    }

    Arborescence res = g.solve(root);
    if (cnt > 1 || !res.exists) write(-1);
    else write(res.weight);
    return 0;
}
