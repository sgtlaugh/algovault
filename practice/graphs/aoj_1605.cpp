/***
 *
 * Weighted Matroid Intersection
 * Maximum weight set of ground elements that is independent in two matroids at once
 *
 * Complexity:
 *   - weighted_matroid_intersection: O(n r^2) oracle queries, r + 1 builds, O(n^2 r^2) time, for an answer of size r
 *   - GraphicMatroid: build O(v + n), queries O(1)
 *   - PartitionMatroid: build O(colors + r), queries O(1)
 *   - XorMatroid: build O(64 n), queries O(1)
 *
 * weighted_matroid_intersection(n, weight, m1, m2): a max weight set of elements 0 to n - 1 independent in both, sorted
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
 *   - weight[x] may be negative, the sum of all |weight[x]| must fit in long long
 *   - Best among the largest sets: add B = 2 sum |weight[x]| + 1 to every weight, n (B + max |weight[x]|) must fit
 *   - Cheapest instead of heaviest: negate every weight
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Augments along the cheapest exchange path, ties broken by fewest elements, which keeps chosen independent in both
template<typename M1, typename M2>
vector<int> weighted_matroid_intersection(int n, const vector<long long>& weight, M1& m1, M2& m2){
    vector<char> in_set(n, 0);
    vector<int> chosen;

    while (true){
        m1.build(chosen), m2.build(chosen);

        vector<vector<int>> next(n);
        for (int y : chosen){
            for (int x = 0; x < n; x++){
                if (in_set[x]) continue;
                if (m1.can_exchange(y, x)) next[y].push_back(x);
                if (m2.can_exchange(y, x)) next[x].push_back(y);
            }
        }

        vector<long long> cost(n, LLONG_MAX);
        vector<int> hops(n, 0), parent(n, -1);
        for (int x = 0; x < n; x++){
            if (!in_set[x] && m1.can_add(x)) cost[x] = -weight[x];
        }
        for (bool changed = true; changed;){
            changed = false;
            for (int u = 0; u < n; u++){
                if (cost[u] == LLONG_MAX) continue;
                for (int v : next[u]){
                    long long c = cost[u] + (in_set[v] ? weight[v] : -weight[v]);
                    if (c < cost[v] || (c == cost[v] && hops[u] + 1 < hops[v])){
                        cost[v] = c, hops[v] = hops[u] + 1, parent[v] = u, changed = true;
                    }
                }
            }
        }

        int sink = -1;
        for (int x = 0; x < n; x++){
            if (in_set[x] || cost[x] == LLONG_MAX || !m2.can_add(x)) continue;
            if (sink == -1 || cost[x] < cost[sink] || (cost[x] == cost[sink] && hops[x] < hops[sink])) sink = x;
        }
        /// The best weight of a k element set is concave in k, so the first path that gains nothing ends the search
        if (sink == -1 || cost[sink] >= 0) break;

        for (int x = sink; x != -1; x = parent[x]) in_set[x] ^= 1;
        chosen.clear();
        for (int x = 0; x < n; x++){
            if (in_set[x]) chosen.push_back(x);
        }
    }

    return chosen;
}

/// BEGIN COPY matroid_oracles from code_library/graphs/matroid_intersection.cpp
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
/// END COPY matroid_oracles

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
 *   - A - not followed by a digit reads as 0 and eats the byte after it, so "3 - 4" reads 3, then 0, then 4
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
 *   - A float, double or long double read gives strtof's, strtod's or strtold's value, 2-3x faster on short decimals
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
        while (isdigit(c)){
            v = v * 10 - (c - '0');

            /// The digits already buffered are read through a local pointer, 15% faster than a read_char call per digit
            const char *p = inbuf + inptr, *end = inbuf + buf_len;
            while (p < end && (unsigned)(*p - '0') < 10) v = v * 10 - (*p++ - '0');
            inptr = p - inbuf;

            /// A buffered non-digit ends the number and is eaten, else the digits may go on after a refill
            if (p < end){
                inptr++;
                break;
            }
            c = read_char();
        }
        x = neg ? v : -v;

        return true;
    }

    /// The C locale's isspace set, false for EOF, and 1.4-1.7x faster for strings and chars than isspace's table lookup
    inline bool is_space(int c){
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (is_space(c)) c = read_char();
        if (c == EOF) return 0;

        for (; !is_space(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return len;
    }

    /// Digits m and 10^|e| both exact in F round once, as parse does: 2-3x faster than parse on short tokens
    template <typename F, F (*parse)(const char *, char **)>
    bool read_floating(F &x, int max_len, unsigned long long max_mantissa, int max_exp){
        int len = read_token();
        if (!len) return false;

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds float and double twice, long double once
        bool own_width = !FLT_EVAL_METHOD || is_same<F, long double>::value;

        /// A longer token skips the digit scan, which would make its parse about 10% slower
        if (len <= max_len && own_width){
            static const F p10[] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L,
                                    1e10L, 1e11L, 1e12L, 1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L,
                                    1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};
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

            if (!*s && digits && digits <= 19 && m <= max_mantissa && e >= -max_exp && e <= max_exp){
                F v = e < 0 ? m / p10[-e] : m * p10[e];
                x = *tmpbuf == '-' ? -v : v;
                return true;
            }
        }

        x = parse(tmpbuf, nullptr);
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        return read_floating<float, strtof>(x, 10, 1 << 24, 10);
    }

    /// 22 characters keep 16-digit scientific tokens like -1.234567890123456e-05 on the fast path
    bool read_one(double &x){
        return read_floating<double, strtod>(x, 22, 1ULL << 53, 22);
    }

    /// Only an x87 64-bit mantissa takes the fast path: a double-double long double does not round a division once
    bool read_one(long double &x){
        return read_floating<long double, strtold>(x, 21, ~0ULL, LDBL_MANT_DIG == 64 ? 27 : -1);
    }

    bool read_one(string &s){
        int c = ' ';

        while (is_space(c)) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !is_space(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (is_space(c));
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

#define MAX 610

int n, m, k, U[MAX], V[MAX], W[MAX], company[MAX], parent[MAX];

int find(int x){
    return parent[x] == x ? x : (parent[x] = find(parent[x]));
}

bool by_cost(int a, int b){
    return W[a] < W[b];
}

vector<int> spanning_forest(int c){
    vector<int> order, res;

    for (int i = 0; i < m; i++){
        if (company[i] == c) order.push_back(i);
    }
    sort(order.begin(), order.end(), by_cost);

    for (int i = 0; i < n; i++) parent[i] = i;
    for (int i : order){
        int x = find(U[i]), y = find(V[i]);
        if (x != y) parent[x] = y, res.push_back(i);
    }

    return res;
}

int main(){
    int i;
    char c;

    while (read(n, m, k) == 3 && (n || m || k)){
        for (i = 0; i < m; i++){
            read(U[i], V[i], W[i], c);
            U[i]--, V[i]--;
            company[i] = (c == 'B');
        }

        vector<int> edges = spanning_forest(0), other = spanning_forest(1);
        edges.insert(edges.end(), other.begin(), other.end());

        int r = edges.size();
        long long total = 0, res = 0;
        GraphicMatroid forest(n);
        vector<int> color(r);
        vector<long long> weight(r);
        for (i = 0; i < r; i++){
            forest.add_edge(U[edges[i]], V[edges[i]]);
            color[i] = company[edges[i]];
            total += W[edges[i]];
        }

        for (i = 0; i < r; i++) weight[i] = (2 * total + 1) - W[edges[i]];
        PartitionMatroid split(color, {k, n - 1 - k});

        vector<int> chosen = weighted_matroid_intersection(r, weight, forest, split);
        for (int x : chosen) res += W[edges[x]];
        if ((int)chosen.size() != (n - 1)) res = -1;
        write(res);
    }

    return 0;
}
