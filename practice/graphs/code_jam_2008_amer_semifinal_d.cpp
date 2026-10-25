/***
 *
 * General Graph Matching
 * Maximum matching in a general graph: Blossom (Edmonds) gives the pairs, TutteMatching (randomized) only the size
 *
 * Complexity:
 *   - Blossom: O(n (n^2 + m)), so O(n^3) without repeated edges, O(n + m) memory
 *   - TutteMatching: O(n^3) even on sparse graphs, O(n^2) memory
 *
 * Blossom b(n): vertices 0 to n - 1, no edges
 * b.add_edge(u, v): undirected edge, u != v, so skip self loops when reading the input
 * b.maximum_matching(): the size of a maximum matching, callable again after adding edges
 * b.mate[v]: after maximum_matching(), the vertex matched to v or -1
 *
 * TutteMatching g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v): undirected edge, u != v
 * g.maximum_matching(): the maximum matching size with probability >= 1 - n / mod, a wrong answer is too small
 *
 * Notes:
 *   - Use Blossom: exact and faster in practice, TutteMatching is a size-only cross-check
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TutteMatching{
    static constexpr unsigned int mod = 1073750017; /// prime, constexpr so % compiles to a multiply

    int n;
    vector<vector<bool>> adj;
    vector<vector<unsigned int>> tutte_matrix;
    mt19937 rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());

    TutteMatching(int n): n(n){
        adj = vector<vector<bool>>(n, vector<bool>(n, 0));
        tutte_matrix = vector<vector<unsigned int>>(n, vector<unsigned int>(n, 0));
    }

    void add_edge(int u, int v){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        adj[u][v] = adj[v][u] = 1;
    }

    int maximum_matching(){
        build_matrix();
        return get_matrix_rank() / 2;
    }

    void build_matrix(){
        for (auto &v : tutte_matrix) fill(v.begin(), v.end(), 0);
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                if (adj[i][j]){
                    unsigned int v = rng() % (mod - 1) + 1;
                    tutte_matrix[i][j] = v, tutte_matrix[j][i] = mod - v;
                }
            }
        }
    }

    int get_matrix_rank(){
        int r = 0;
        vector<int> nonzero;

        for (int j = 0; j < n; j++){
            int k = r;
            while (k < n && !tutte_matrix[k][j]) k++;
            if (k == n) continue;

            swap(tutte_matrix[k], tutte_matrix[r]);
            auto& pivot = tutte_matrix[r];
            unsigned long long inv = expo(pivot[j], mod - 2);

            nonzero.clear();
            for (int v = j + 1; v < n; v++){
                if (pivot[v]){
                    pivot[v] = inv * pivot[v] % mod;
                    nonzero.push_back(v);
                }
            }

            for (int u = r + 1; u < n; u++){
                auto& row = tutte_matrix[u];
                if (!row[j]) continue;
                unsigned long long f = mod - row[j];
                for (int v : nonzero) row[v] = (row[v] + f * pivot[v]) % mod;
            }
            r++;
        }

        return r;
    }

    unsigned int expo(unsigned long long x, unsigned int n){
        unsigned long long res = 1;

        while (n){
            if (n & 1) res = res * x % mod;
            x = x * x % mod;
            n >>= 1;
        }

        return res;
    }
};

struct Blossom{
    int n;
    vector<vector<int>> adj;
    vector<int> mate, parent, base, seen_stamp;
    vector<bool> in_tree, in_blossom;
    int stamp = 0;

    Blossom(int n) : n(n), adj(n), mate(n, -1), parent(n), base(n), seen_stamp(n, 0), in_tree(n), in_blossom(n) {}

    void add_edge(int u, int v){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    int maximum_matching(){
        fill(mate.begin(), mate.end(), -1);
        for (int u = 0; u < n; u++){
            for (int v : adj[u]){
                if (mate[u] == -1 && mate[v] == -1) mate[u] = v, mate[v] = u;
            }
        }

        for (int root = 0; root < n; root++){
            if (mate[root] != -1) continue;
            for (int v = find_augmenting_path(root); v != -1;){
                int next = mate[parent[v]];
                mate[v] = parent[v], mate[parent[v]] = v;
                v = next;
            }
        }

        int size = 0;
        for (int v = 0; v < n; v++) size += mate[v] > v;
        return size;
    }

    /// Returns the free vertex ending an augmenting path from root, parent[] leads back along it, or -1 if none exists
    int find_augmenting_path(int root){
        fill(parent.begin(), parent.end(), -1);
        fill(in_tree.begin(), in_tree.end(), false);
        iota(base.begin(), base.end(), 0);
        /// A search contracts at most n / 2 blossoms, so resetting here keeps stamp far from overflow
        stamp = 0;
        fill(seen_stamp.begin(), seen_stamp.end(), 0);

        vector<int> queue = {root};
        in_tree[root] = true;
        for (size_t i = 0; i < queue.size(); i++){
            int v = queue[i];
            for (int to : adj[v]){
                if (base[v] == base[to] || mate[v] == to) continue;

                /// to is an even (outer) vertex of the tree, the edge closes an odd cycle
                if (to == root || (mate[to] != -1 && parent[mate[to]] != -1)){
                    int lca = lowest_common_base(v, to);
                    fill(in_blossom.begin(), in_blossom.end(), false);
                    mark_path(v, lca, to);
                    mark_path(to, lca, v);
                    for (int u = 0; u < n; u++){
                        if (!in_blossom[base[u]]) continue;
                        base[u] = lca;
                        if (!in_tree[u]) in_tree[u] = true, queue.push_back(u);
                    }
                }
                else if (parent[to] == -1){
                    parent[to] = v;
                    if (mate[to] == -1) return to;
                    in_tree[mate[to]] = true;
                    queue.push_back(mate[to]);
                }
            }
        }

        return -1;
    }

    int lowest_common_base(int u, int v){
        stamp++;
        while (true){
            u = base[u];
            seen_stamp[u] = stamp;
            if (mate[u] == -1) break;
            u = parent[mate[u]];
        }

        while (true){
            v = base[v];
            if (seen_stamp[v] == stamp) return v;
            v = parent[mate[v]];
        }
    }

    /// Walks from v up to base lca, flagging the bases it crosses and pointing odd vertices back across the cycle
    void mark_path(int v, int lca, int child){
        while (base[v] != lca){
            in_blossom[base[v]] = in_blossom[base[mate[v]]] = true;
            parent[v] = child;
            child = mate[v];
            v = parent[mate[v]];
        }
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

int r, c;
string grid[16];

int matching(bool with_king){
    int i, j, d, e, x, y;
    Blossom b(r * c);

    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            if (grid[i][j] == '#' || (!with_king && grid[i][j] == 'K')) continue;
            for (d = -1; d <= 1; d++){
                for (e = -1; e <= 1; e++){
                    x = i + d, y = j + e;
                    if ((d == 0 && e == 0) || x < 0 || x >= r || y < 0 || y >= c) continue;
                    if (grid[x][y] == '#' || (!with_king && grid[x][y] == 'K')) continue;
                    if ((x * c + y) > (i * c + j)) b.add_edge(i * c + j, x * c + y);
                }
            }
        }
    }

    return b.maximum_matching();
}

int main(){
    int t, T = 0, i;
    char buf[32];

    read(t);
    while (t--){
        read(r, c);
        for (i = 0; i < r; i++) read(grid[i]);

        snprintf(buf, sizeof(buf), "Case #%d: %c", ++T, matching(true) > matching(false) ? 'A' : 'B');
        write(string(buf));
    }
    return 0;
}
