/***
 *
 * Matroid Intersection (unweighted)
 * Largest set of ground elements that is independent in two matroids at once
 *
 * Complexity:
 *   - matroid_intersection: O(n r^2 + n) oracle queries and r + 1 builds, for n elements and an answer of size r
 *   - GraphicMatroid: build O(v + n), queries O(1)
 *   - PartitionMatroid: build O(colors + r), queries O(1)
 *   - XorMatroid: build O(64 n), queries O(1)
 *
 * matroid_intersection(n, m1, m2): the largest set of elements 0 to n - 1 independent in both, sorted ascending
 * m1 and m2 are oracles, any structs with these three calls:
 *     void build(const vector<int>& chosen): chosen is the current common independent set, called before each round
 *     bool can_add(int x): true if chosen + x is independent, x not in chosen
 *     bool can_exchange(int y, int x): true if chosen - y + x is independent, y in chosen, x not in chosen
 *
 * GraphicMatroid g(v): independent sets are forests on vertices 0 to v - 1, self-loops are never independent
 * g.add_edge(a, b): adds the edge a - b as element g.edges.size()
 * PartitionMatroid p(color, cap): element x has color[x], an independent set has at most cap[c] elements of color c
 * XorMatroid x(vec): independent sets of 64-bit vectors are linearly independent over GF(2), 0 never is
 *
 * Notes:
 *   - A rainbow spanning tree is GraphicMatroid intersected with a PartitionMatroid whose caps are all 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename M1, typename M2>
vector<int> matroid_intersection(int n, M1& m1, M2& m2){
    const int SOURCE = -1, UNSEEN = -2;
    vector<char> in_set(n, 0);
    vector<int> chosen;

    while (true){
        m1.build(chosen), m2.build(chosen);

        vector<int> parent(n, UNSEEN), bfs;
        for (int x = 0; x < n; x++){
            if (!in_set[x] && m1.can_add(x)) parent[x] = SOURCE, bfs.push_back(x);
        }

        int sink = -1;
        for (int i = 0; i < (int)bfs.size() && sink == -1; i++){
            int u = bfs[i];
            if (in_set[u]){
                for (int x = 0; x < n; x++){
                    if (!in_set[x] && parent[x] == UNSEEN && m1.can_exchange(u, x)) parent[x] = u, bfs.push_back(x);
                }
            }
            else if (m2.can_add(u)) sink = u;
            else{
                for (int y : chosen){
                    if (parent[y] == UNSEEN && m2.can_exchange(y, u)) parent[y] = u, bfs.push_back(y);
                }
            }
        }
        if (sink == -1) break;

        for (int x = sink; x != SOURCE; x = parent[x]) in_set[x] ^= 1;
        chosen.clear();
        for (int x = 0; x < n; x++){
            if (in_set[x]) chosen.push_back(x);
        }
    }

    return chosen;
}

/// BEGIN SHARED matroid_oracles
/// chosen - y + x is a forest iff x joins two trees of chosen, or y lies on the tree path between x's endpoints
struct GraphicMatroid{
    int v;
    vector<array<int, 2>> edges;
    vector<vector<int>> adj;
    vector<int> comp, tin, tout, child, it;

    GraphicMatroid(int v) : v(v) {}

    void add_edge(int a, int b){
        edges.push_back({a, b});
    }

    void build(const vector<int>& chosen){
        adj.assign(v, {});
        for (int e : chosen){
            adj[edges[e][0]].push_back(e);
            adj[edges[e][1]].push_back(e);
        }

        comp.assign(v, -1), tin.assign(v, 0), tout.assign(v, 0), it.assign(v, 0);
        child.assign(edges.size(), -1);
        int timer = 0;
        for (int root = 0; root < v; root++){
            if (comp[root] != -1) continue;
            comp[root] = root, tin[root] = timer++;
            vector<int> stack = {root};
            while (!stack.empty()){
                int a = stack.back();
                if (it[a] == (int)adj[a].size()){
                    tout[a] = timer;
                    stack.pop_back();
                    continue;
                }

                int e = adj[a][it[a]++], b = edges[e][0] ^ edges[e][1] ^ a;
                if (comp[b] != -1) continue;
                comp[b] = root, tin[b] = timer++, child[e] = b;
                stack.push_back(b);
            }
        }
    }

    bool can_add(int x) const{
        return comp[edges[x][0]] != comp[edges[x][1]];
    }

    bool can_exchange(int y, int x) const{
        if (can_add(x)) return true;
        int c = child[y];
        return inside(c, edges[x][0]) != inside(c, edges[x][1]);
    }

    bool inside(int c, int a) const{
        return tin[c] <= tin[a] && tin[a] < tout[c];
    }
};

struct PartitionMatroid{
    vector<int> color, cap, used;

    PartitionMatroid(const vector<int>& color, const vector<int>& cap) : color(color), cap(cap) {}

    void build(const vector<int>& chosen){
        used.assign(cap.size(), 0);
        for (int x : chosen) used[color[x]]++;
    }

    bool can_add(int x) const{
        return used[color[x]] < cap[color[x]];
    }

    bool can_exchange(int y, int x) const{
        return color[x] == color[y] || can_add(x);
    }
};

/// Writes each element over the basis of chosen: chosen - y + x is independent iff x is outside span(chosen) or uses y
struct XorMatroid{
    vector<unsigned long long> vec, uses;
    vector<char> spanned;
    vector<int> pos;

    XorMatroid(const vector<unsigned long long>& vec) : vec(vec) {}

    void build(const vector<int>& chosen){
        int n = vec.size();
        unsigned long long basis[64] = {}, basis_uses[64] = {};
        pos.assign(n, -1);
        for (int i = 0; i < (int)chosen.size(); i++){
            unsigned long long a = vec[chosen[i]], mask = 1ULL << i;
            pos[chosen[i]] = i;
            for (int b = 63; b >= 0; b--){
                if (!(a >> b & 1)) continue;
                if (!basis[b]){
                    basis[b] = a, basis_uses[b] = mask;
                    break;
                }
                a ^= basis[b], mask ^= basis_uses[b];
            }
        }

        uses.assign(n, 0), spanned.assign(n, 0);
        for (int x = 0; x < n; x++){
            unsigned long long a = vec[x], mask = 0;
            for (int b = 63; b >= 0; b--){
                if ((a >> b & 1) && basis[b]) a ^= basis[b], mask ^= basis_uses[b];
            }
            spanned[x] = a == 0, uses[x] = mask;
        }
    }

    bool can_add(int x) const{
        return !spanned[x];
    }

    bool can_exchange(int y, int x) const{
        return !spanned[x] || (uses[x] >> pos[y] & 1);
    }
};
/// END SHARED matroid_oracles

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

unsigned long long basis[64];

bool insert(unsigned long long x){
    for (int b = 63; b >= 0; b--){
        if (!((x >> b) & 1)) continue;
        if (!basis[b]){
            basis[b] = x;
            return true;
        }
        x ^= basis[b];
    }

    return false;
}

unsigned long long reduce(unsigned long long x){
    for (int b = 63; b >= 0; b--){
        if (((x >> b) & 1) && basis[b]) x ^= basis[b];
    }

    return x;
}

int main(){
    int n, m, i, j, k;
    unsigned long long x;

    read(n);
    bool flag = true;
    for (i = 0; i < n; i++){
        read(x);
        if (!insert(x)) flag = false;
    }

    read(m);
    vector<unsigned long long> heaps, vec;
    vector<int> box, cap(m, 1);
    for (i = 0; i < m; i++){
        read(k);
        for (j = 0; j < k; j++){
            read(x);
            heaps.push_back(x), vec.push_back(reduce(x)), box.push_back(i);
        }
    }

    if (!flag){
        write(-1);
        return 0;
    }

    XorMatroid m1(vec);
    PartitionMatroid m2(box, cap);
    vector<int> res = matroid_intersection(heaps.size(), m1, m2);
    if ((int)res.size() < m){
        write(-1);
        return 0;
    }

    vector<unsigned long long> pick(m);
    for (int e: res) pick[box[e]] = heaps[e];
    for (i = 0; i < m; i++) write(pick[i]);
    return 0;
}
