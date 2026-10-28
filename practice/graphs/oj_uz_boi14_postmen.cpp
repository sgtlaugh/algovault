/***
 *
 * Euler Path (Hierholzer)
 * Finds a walk that uses every edge exactly once, in a directed or undirected multigraph
 *
 * Complexity: O(n + m)
 *
 * EulerDirected g(n) or EulerUndirected g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v): adds an edge and returns its id, ids are 0, 1, 2, ... in call order
 * g.solve(start): the m + 1 vertices of an Euler path from start, empty if none exists, start is -1 or a vertex
 * g.walk_edges: after a successful solve, the m edge ids in order, edge walk_edges[i] joins path[i] to path[i + 1]
 *
 * Notes:
 *   - start = -1 picks the first vertex with out - in = 1 (directed) or odd degree, else the first with an edge
 *   - The path is a cycle iff path.front() == path.back()
 *   - With no edges the path is the single vertex start (0 for start = -1), and it is empty when n = 0
 *   - Multi-edges, self loops and isolated vertices are allowed: only the vertices that have edges need to be connected
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<bool directed>
struct EulerPath{
    int n, m = 0;
    vector<vector<pair<int, int>>> adj;
    vector<int> balance, walk_edges; /// balance: out - in when directed, the degree when undirected

    EulerPath(int n) : n(n), adj(n), balance(n, 0) {}

    int add_edge(int u, int v){
        adj[u].push_back({v, m});
        if (directed) balance[u]++, balance[v]--;
        else adj[v].push_back({u, m}), balance[u]++, balance[v]++;
        return m++;
    }

    vector<int> solve(int start = -1){
        walk_edges.clear();
        if (n == 0) return {};

        int unbalanced = 0, source = -1;
        for (int u = 0; u < n; u++){
            int b = excess(u);
            if (b == 0) continue;
            if (b != 1 && b != -1) return {};
            unbalanced++;
            if (b == 1 && source == -1) source = u;
        }
        if (unbalanced > 2) return {};

        if (start == -1){
            start = unbalanced ? source : 0;
            for (int u = 0; u < n && !unbalanced; u++){
                if (!adj[u].empty()){
                    start = u;
                    break;
                }
            }
        }
        if (unbalanced ? excess(start) != 1 : m > 0 && adj[start].empty()) return {};

        vector<int> it(n, 0), path;
        vector<char> used(m, 0);
        vector<pair<int, int>> stack = {{start, -1}};
        while (!stack.empty()){
            auto [u, e] = stack.back();
            if (it[u] == (int)adj[u].size()){
                path.push_back(u);
                if (e != -1) walk_edges.push_back(e);
                stack.pop_back();
                continue;
            }
            auto [v, id] = adj[u][it[u]++];
            if (!used[id]) used[id] = 1, stack.push_back({v, id});
        }

        if ((int)path.size() != m + 1){
            walk_edges.clear();
            return {};
        }
        reverse(path.begin(), path.end());
        reverse(walk_edges.begin(), walk_edges.end());
        return path;
    }

    int excess(int u) const{
        return directed ? balance[u] : balance[u] % 2;
    }
};

using EulerDirected = EulerPath<true>;
using EulerUndirected = EulerPath<false>;

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

#define MAX 500010

int ea[MAX], eb[MAX], deg[MAX];
bool in_stack[MAX];

int main(){
    int n, m, i;

    read(n, m);
    for (i = 0; i < m; i++){
        read(ea[i], eb[i]);
        deg[--ea[i]]++, deg[--eb[i]]++;
    }

    EulerUndirected g(n);
    for (i = 0; i < n; i++) g.adj[i].reserve(deg[i]);
    for (i = 0; i < m; i++) g.add_edge(ea[i], eb[i]);

    vector<int> path = g.solve(0), st;
    for (int v: path){
        if (in_stack[v]){
            vector<int> tour = {v + 1};
            while (st.back() != v){
                tour.push_back(st.back() + 1);
                in_stack[st.back()] = false;
                st.pop_back();
            }
            write(tour);
        }
        else st.push_back(v), in_stack[v] = true;
    }

    return 0;
}
