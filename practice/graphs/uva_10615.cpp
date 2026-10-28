/***
 *
 * Edge Coloring
 * Colors the edges so that edges sharing a vertex get different colors
 *
 * Complexity:
 *   - BipartiteEdgeColoring: O(n + m n) expected time
 *   - VizingEdgeColoring: O(n m) time, O(n D) memory for maximum degree D
 *
 * BipartiteEdgeColoring g(n_left, n_right): left vertices 0..n_left - 1, right vertices 0..n_right - 1
 * g.add_edge(u, v): edge from left u to right v, parallel edges allowed
 * g.solve(): color of every edge in insertion order, in 0..D - 1, optimal by Konig's theorem
 *
 * VizingEdgeColoring g(n): vertices 0..n - 1, the graph must be simple (no loops, no parallel edges)
 * g.add_edge(u, v): undirected edge
 * g.solve(): color of every edge in insertion order, in 0..D (Misra-Gries), as a D coloring is NP-hard in general
 *
 * Notes:
 *   - Both solve() calls can be repeated after adding more edges
 *   - VizingEdgeColoring allocates n (D + 1) ints, so a star with n = 1e5 needs 40 GB
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct BipartiteEdgeColoring{
    int n_left, n_right;
    vector<array<int, 2>> edges;

    BipartiteEdgeColoring(int n_left, int n_right) : n_left(n_left), n_right(n_right) {}

    void add_edge(int u, int v){
        edges.push_back({u, v});
    }

    vector<int> solve(){
        int m = edges.size(), k = 0;
        vector<int> deg[2] = {vector<int>(n_left), vector<int>(n_right)};
        for (auto& e : edges){
            k = max({k, ++deg[0][e[0]], ++deg[1][e[1]]});
        }

        /// A free color at x always exists below deg(x): candidates[s][x] holds every free one, stale ones pop lazily
        vector<vector<int>> candidates[2];
        for (int s = 0; s < 2; s++){
            for (int d : deg[s]){
                candidates[s].emplace_back(d);
                iota(candidates[s].back().rbegin(), candidates[s].back().rend(), 0);
            }
        }

        /// (side, vertex, color) -> edge or -1, hashed instead of an O(n D) table, 8m slots for at most 4m keys ever
        int bits = 1;
        while ((1LL << bits) < 8LL * m) bits++;
        vector<pair<long long, int>> table(1LL << bits, {-1, -1});
        vector<int> color(m, -1), path;
        auto slot = [&](int s, int x, int c) -> int& {
            long long key = (s ? (long long)n_left + x : x) * k + c;
            size_t i = (unsigned long long)key * 0x9E3779B97F4A7C15ULL >> (64 - bits);
            while (table[i].first != key && table[i].first != -1) i = (i + 1) & (table.size() - 1);
            if (table[i].first == -1) table[i] = {key, -1};
            return table[i].second;
        };
        auto first_free = [&](int s, int x){
            auto& list = candidates[s][x];
            while (slot(s, x, list.back()) != -1) list.pop_back();
            return list.back();
        };

        for (int e = 0; e < m; e++){
            int u = edges[e][0], v = edges[e][1], a = first_free(0, u), b = first_free(1, v);
            if (slot(1, v, a) != -1){
                /// The a/b path from v enters left vertices by a edges, so it never reaches u where a is free
                path.clear();
                int side = 1, x = v, c = a;
                for (; slot(side, x, c) != -1; side ^= 1, c ^= a ^ b){
                    int f = slot(side, x, c);
                    path.push_back(f);
                    x = edges[f][side ^ 1];
                }
                /// The path's far end loses color c ^ a ^ b, the only color the inversion frees apart from a at v
                candidates[side][x].push_back(c ^ a ^ b);

                for (int f : path){
                    for (int s = 0; s < 2; s++) slot(s, edges[f][s], color[f]) = -1;
                }
                for (int f : path){
                    color[f] ^= a ^ b;
                    for (int s = 0; s < 2; s++) slot(s, edges[f][s], color[f]) = f;
                }
            }
            color[e] = a;
            slot(0, u, a) = slot(1, v, a) = e;
        }

        return color;
    }
};

struct VizingEdgeColoring{
    int n;
    vector<array<int, 2>> edges;

    VizingEdgeColoring(int n) : n(n) {}

    void add_edge(int u, int v){
        assert(u != v);
        edges.push_back({u, v});
    }

    vector<int> solve(){
        int m = edges.size(), k = 1;
        vector<int> deg(n);
        for (auto& e : edges){
            k = max({k, ++deg[e[0]] + 1, ++deg[e[1]] + 1});
        }

        /// adj[x * k + c]: neighbor of x through color c or -1, flat since hashing it ran 14x slower on dense graphs
        vector<int> adj((size_t)n * k, -1), free_color(n), fan(k + 1), fan_color(k + 1), seen;
        auto at = [&](int x, int c) -> int& { return adj[(long long)x * k + c]; };

        for (auto& e : edges){
            int u = e[0], v = e[1], c = free_color[u], d, len = 0;
            fan[0] = v;
            seen.assign(k, 0);
            while (d = free_color[v], !seen[d] && (v = at(u, d)) != -1){
                seen[d] = ++len;
                fan_color[len] = d, fan[len] = v;
            }

            /// Inverts the c/d path starting at u, which recolors u's d edge (fan[seen[d]], if any) to c
            fan_color[seen[d]] = c;
            int end = u;
            for (int x = u, cd = d; x != -1; cd ^= c ^ d, x = at(x, cd)){
                swap(at(x, cd), at(x, cd ^ c ^ d));
                end = x;
            }

            /// Rotates the fan until its current tip has d free, then colors the tip's edge d
            int i = 0;
            while (at(fan[i], d) != -1){
                int left = fan[i], right = fan[++i], f = fan_color[i];
                at(u, f) = left, at(left, f) = u;
                at(right, f) = -1;
                free_color[right] = f;
            }
            at(u, d) = fan[i], at(fan[i], d) = u;

            for (int y : {fan[0], u, end}){
                free_color[y] = 0;
                while (at(y, free_color[y]) != -1) free_color[y]++;
            }
        }

        vector<int> color(m, 0);
        for (int i = 0; i < m; i++){
            while (at(edges[i][0], color[i]) != edges[i][1]) color[i]++;
        }

        return color;
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

#define MAX 110

int n, color[MAX][MAX];
string board[MAX];

int main(){
    int t, i, j, k, res;

    while (read(t)){
        while (t--){
            read(n);
            for (i = 0; i < n; i++) read(board[i]);

            BipartiteEdgeColoring g(n, n);
            for (i = 0; i < n; i++){
                for (j = 0; j < n; j++){
                    if (board[i][j] == '*') g.add_edge(i, j);
                }
            }
            vector<int> colors = g.solve();

            k = 0, res = 0;
            for (i = 0; i < n; i++){
                for (j = 0; j < n; j++){
                    color[i][j] = (board[i][j] == '*') ? colors[k++] + 1 : 0;
                    res = max(res, color[i][j]);
                }
            }

            write(res);
            for (i = 0; i < n; i++) write(vector<int>(color[i], color[i] + n));
        }
    }

    return 0;
}
