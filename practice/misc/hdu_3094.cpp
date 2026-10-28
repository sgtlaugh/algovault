/***
 *
 * Hackenbush
 * Green Hackenbush Grundy value of a rooted graph and red-blue Hackenbush value of a rooted tree
 *
 * Complexity:
 *   - GreenHackenbush grundy: O(n + m) per call
 *   - RedBlueHackenbush value: O(n log^2 n + m) per call
 *
 * GreenHackenbush g(n): impartial game on nodes 0 to n - 1, either player deletes any edge
 * g.add_edge(u, v): parallel edges and self loops (u == v) allowed
 * g.grundy(root): the Grundy value, the first player wins iff it is nonzero
 *
 * RedBlueHackenbush t(n): partizan game on nodes 0 to n - 1, Left deletes only blue edges, Right only red
 * t.add_edge(u, v, blue): adds a blue edge if blue, else a red one
 * t.value(root): the exact value as a Dyadic, the root's component must be a tree
 *
 * Dyadic: whole + sum of 2^-b over b in bits, whole = floor, bits strictly increasing, all >= 1
 * d.sign(): > 0 if Left wins whoever starts, < 0 if Right wins whoever starts, 0 if the second player wins
 *
 * Notes:
 *   - A move deletes one edge and every edge it cuts off from the root, a player who cannot move loses
 *   - Several ground vertices: map them all to one root vertex before add_edge
 *   - Edges outside the root's component take no part in the game
 *   - In RedBlueHackenbush parallel edges and self loops are cycles
 *   - Denominators reach 2^(n - 2) on a long stalk, so the value is kept as bits, not as a double
 *   - Fusion principle: the vertices of a 2-edge-connected part fuse into one, its edges become loops (1 each)
 *   - Colon principle: a subtree of value x on top of a bridge is worth x + 1, siblings combine by xor
 *   - A blue edge over a subtree of value x is worth x + 1 for x >= 0, else (1 + frac(x)) / 2^(-floor(x))
 *   - A red edge over x is worth the negation of the blue value of -x, and sibling values add
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GreenHackenbush{
    struct Frame{
        int u, up, next, end, dis, low, value;
    };

    int n;
    vector<pair<int, int>> edges;

    GreenHackenbush(int n) : n(n) {}

    void add_edge(int u, int v){
        edges.push_back({u, v});
    }

    int grundy(int root){
        vector<int> degree(n, 0), neighbors(n, 0), value(n, 0); /// neighbors xors the ids, the only one at degree 1
        for (auto e : edges){
            int u = e.first, v = e.second;
            if (u == v){
                value[u] ^= 1;
                continue;
            }
            degree[u]++, neighbors[u] ^= v;
            degree[v]++, neighbors[v] ^= u;
        }

        /// Pendant trees fold into the vertex they hang from by the colon principle, with no adjacency lists built
        vector<int> leaves;
        for (int u = 0; u < n; u++){
            if (u != root && degree[u] == 1) leaves.push_back(u);
        }
        for (int i = 0; i < (int)leaves.size(); i++){
            int u = leaves[i], p = neighbors[u];
            if (degree[u] != 1) continue; /// its last neighbor was folded into it, in a component without the root

            degree[u] = 0;
            degree[p]--, neighbors[p] ^= u, value[p] ^= value[u] + 1;
            if (p != root && degree[p] == 1) leaves.push_back(p);
        }
        if (degree[root] == 0) return value[root];

        /// What is left is cycles and paths between them: Tarjan bridges on the edges with both ends unfolded
        vector<int> start(n + 1, 0), adj, dis(n, -1);
        for (auto e : edges){
            int u = e.first, v = e.second;
            if (degree[u] && degree[v]) start[u]++, start[v]++;
        }
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        adj.resize(start[n]);
        for (auto e : edges){
            int u = e.first, v = e.second;
            if (degree[u] && degree[v]) adj[--start[u]] = v, adj[--start[v]] = u;
        }

        int timer = 0;
        vector<Frame> stack = {{root, -1, start[root], start[root + 1], timer, timer, value[root]}};
        dis[root] = timer++;

        while (true){
            Frame& f = stack.back();
            if (f.next < f.end){
                int v = adj[f.next++];
                /// One edge to the parent is the tree edge, any parallel copy of it closes a cycle
                if (v == f.up) f.up = -1;
                else if (dis[v] == -1){
                    dis[v] = timer;
                    stack.push_back({v, f.u, start[v], start[v + 1], timer, timer, value[v]});
                    timer++;
                }
                else if (dis[v] < f.dis){
                    f.low = min(f.low, dis[v]);
                    f.value ^= 1; /// a cycle edge fuses into a loop worth 1, counted once, at its deeper endpoint
                }
                continue;
            }

            if (stack.size() == 1) return f.value;
            Frame& p = stack[stack.size() - 2];
            p.low = min(p.low, f.low);
            p.value ^= f.low > p.dis ? f.value + 1 : f.value ^ 1;
            stack.pop_back();
        }
    }
};

struct Dyadic{
    long long whole = 0;
    vector<int> bits;

    int sign() const{
        if (whole != 0) return whole > 0 ? 1 : -1;
        return !bits.empty();
    }

    bool operator==(const Dyadic& other) const{
        return whole == other.whole && bits == other.bits;
    }
};

struct RedBlueHackenbush{
    struct Number{
        long long whole = 0;
        int shift = 0;
        set<int> bits; /// actual exponent = stored + shift, so halving every bit at once is O(1)

        void add_bit(int e){
            while (e > 0 && bits.erase(e - shift)) e--;
            if (e == 0) whole++;
            else bits.insert(e - shift);
        }

        void absorb(const Number& other){
            whole += other.whole;
            for (int b : other.bits) add_bit(b + other.shift);
        }

        void below_blue(){
            if (whole >= 0){
                whole++;
                return;
            }

            int k = -whole;
            shift += k, whole = 0;
            bits.insert(k - shift);
        }

        /// (frac - 2) / 2^(w + 1) = -1 + (1 - 2^-w) + frac / 2^(w + 1): w new bits, paid by the w integer part edges
        void below_red(){
            if (whole < 0 || (whole == 0 && bits.empty())){
                whole--;
                return;
            }

            int w = whole;
            shift += w + 1, whole = -1;
            for (int e = 1; e <= w; e++) bits.insert(e - shift);
        }
    };

    int n;
    vector<array<int, 3>> edges; /// (u, v, blue)

    RedBlueHackenbush(int n) : n(n) {}

    void add_edge(int u, int v, bool blue){
        edges.push_back({u, v, blue});
    }

    Dyadic value(int root){
        vector<int> degree(n, 0), incident(n, 0); /// incident xors the edge ids, the only one at degree 1
        for (int id = 0; id < (int)edges.size(); id++){
            int u = edges[id][0], v = edges[id][1], blue = edges[id][2];
            degree[u]++, incident[u] ^= id;
            degree[v]++, incident[v] ^= id;
        }

        /// Leaves fold into their parent until only the root is left, so no adjacency lists are built
        vector<int> leaves;
        for (int u = 0; u < n; u++){
            if (u != root && degree[u] == 1) leaves.push_back(u);
        }

        vector<Number> num(n);
        for (int i = 0; i < (int)leaves.size(); i++){
            int u = leaves[i];
            if (degree[u] != 1) continue; /// its last neighbor was folded into it, in a component without the root

            int a = edges[incident[u]][0], b = edges[incident[u]][1], blue = edges[incident[u]][2];
            int p = a ^ b ^ u;
            degree[u] = 0;
            degree[p]--, incident[p] ^= incident[u];
            if (p != root && degree[p] == 1) leaves.push_back(p);

            blue ? num[u].below_blue() : num[u].below_red();
            if (num[p].bits.size() < num[u].bits.size()) swap(num[p], num[u]);
            num[p].absorb(num[u]);
            num[u] = Number();
        }
        assert(degree[root] == 0); /// a cycle in the root's component never folds away, so the root keeps an edge

        Dyadic result;
        result.whole = num[root].whole;
        for (int b : num[root].bits) result.bits.push_back(b + num[root].shift);
        return result;
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
            GreenHackenbush g(n);
            for (i = 1; i < n; i++){
                read(u, v);
                g.add_edge(u - 1, v - 1);
            }

            if (g.grundy(0)) write("Alice");
            else write("Bob");
        }
    }

    return 0;
}
