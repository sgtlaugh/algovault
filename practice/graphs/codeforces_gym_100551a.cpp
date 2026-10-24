/***
 *
 * Offline Dynamic Connectivity
 * Answers connectivity and component count queries on an undirected graph under edge insertions and deletions
 *
 * Complexity: O(n + (m + q) log(m + q) log n) time, O(n + q + m log q) memory, for m edge operations and q queries
 *
 * DynamicConnectivity dc(n): an empty graph on vertices 0 to n - 1, parallel edges and self loops allowed
 * dc.add_edge(u, v): adds a copy of the edge (u, v)
 * dc.query_components(): records a query for the number of connected components, returns its index
 * dc.query_connected(u, v): records a query for 1 if u and v are connected, else 0, returns its index
 * dc.remove_edge(u, v): deletes one live copy of (u, v) in either orientation, the copy must exist
 * dc.solve(): the answers in query order, indexed by what the query calls returned
 *
 * Notes:
 *   - solve() keeps the operations, so more can be added and solve() called again
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY rollback_dsu from code_library/data_structures/rollback_dsu.cpp
struct RollbackDSU{
    vector<int> sz, parent, history;

    RollbackDSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        history.push_back(a);
        return true;
    }

    void rollback(int t){
        assert(0 <= t && t <= time());
        while ((int)history.size() > t){
            int a = history.back();
            history.pop_back();
            sz[parent[a]] -= sz[a], parent[a] = a;
        }
    }

    int find(int i){
        while (i != parent[i]) i = parent[i];
        return i;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }

    int time(){
        return history.size();
    }
};
/// END COPY rollback_dsu

struct DynamicConnectivity{
    int n;
    vector<array<int, 2>> queries; /// (u, v), or (-1, -1) for a component count
    vector<array<int, 4>> intervals; /// edge (u, v) alive for queries [l, r)
    map<pair<int, int>, vector<int>> open; /// first query index of each live copy of an edge
    vector<vector<pair<int, int>>> tree;
    RollbackDSU dsu;

    DynamicConnectivity(int n): n(n), dsu(n) {}

    void add_edge(int u, int v){
        if (u > v) swap(u, v);
        open[{u, v}].push_back(queries.size());
    }

    int query_components(){
        queries.push_back({-1, -1});
        return queries.size() - 1;
    }

    int query_connected(int u, int v){
        queries.push_back({u, v});
        return queries.size() - 1;
    }

    void remove_edge(int u, int v){
        if (u > v) swap(u, v);
        auto it = open.find({u, v});
        assert(it != open.end());

        int l = it->second.back();
        it->second.pop_back();
        if (it->second.empty()) open.erase(it);
        if (l < (int)queries.size()) intervals.push_back({u, v, l, (int)queries.size()});
    }

    vector<int> solve(){
        int q = queries.size();
        vector<int> res(q);
        if (q == 0) return res;

        tree.assign(4 * q, {});
        for (auto [u, v, l, r]: intervals) insert(1, 0, q - 1, l, r - 1, {u, v});
        for (auto& [edge, starts]: open){
            for (int l: starts){
                if (l < q) insert(1, 0, q - 1, l, q - 1, edge);
            }
        }

        dfs(1, 0, q - 1, res);

        tree.clear();
        return res;
    }

    void insert(int node, int lo, int hi, int l, int r, pair<int, int> edge){
        if (r < lo || hi < l) return;
        if (l <= lo && hi <= r){
            tree[node].push_back(edge);
            return;
        }

        int mid = (lo + hi) / 2;
        insert(2 * node, lo, mid, l, r, edge);
        insert(2 * node + 1, mid + 1, hi, l, r, edge);
    }

    void dfs(int node, int lo, int hi, vector<int>& res){
        int saved = dsu.time();
        for (auto [u, v]: tree[node]) dsu.connect(u, v);

        if (lo == hi){
            auto [u, v] = queries[lo];
            res[lo] = u == -1 ? n - dsu.time() : dsu.same(u, v);
        }
        else{
            int mid = (lo + hi) / 2;
            dfs(2 * node, lo, mid, res);
            dfs(2 * node + 1, mid + 1, hi, res);
        }

        dsu.rollback(saved);
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

int main(){
    int n, k, u, v;
    char c;

    if (fopen("connect.in", "r")){
        freopen("connect.in", "r", stdin);
        freopen("connect.out", "w", stdout);
    }

    read(n, k);
    DynamicConnectivity dc(n);
    while (k--){
        read(c);
        if (c == '?'){
            dc.query_components();
            continue;
        }

        read(u, v);
        if (c == '+') dc.add_edge(u - 1, v - 1);
        else dc.remove_edge(u - 1, v - 1);
    }

    for (int x: dc.solve()) write(x);
    return 0;
}
