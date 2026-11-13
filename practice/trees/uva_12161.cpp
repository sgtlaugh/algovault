/***
 *
 * Centroid Decomposition
 * Recursively cuts a weighted forest at centroids, so the path u-v passes through lca(u, v) in the centroid tree
 *
 * Complexity: O(n log n) to build, plus the visitor's own work on O(n log n) delivered (node, distance) pairs
 *
 * CentroidDecomposition cd(n): an empty forest on nodes 0..n-1
 * cd.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * cd.build(visit): calls visit(c, branches) once per centroid c, parents before children, with c already cut out
 * cd.build(): fills only parent and depth
 * count_paths_at_most(cd, k): the number of unordered pairs u != v with dist(u, v) <= k, in O(n log^2 n)
 * cd.parent[c]: parent of c in the centroid tree, -1 for the root centroid of each tree of the forest
 * cd.depth[c]: depth of c in the centroid tree, at most floor(log2(n))
 *
 * Notes:
 *   - branches[i] lists (node, distance to c) for the piece behind the i-th uncut neighbor of c, one child's component
 *   - branches has type CentroidDecomposition::Branches
 *   - build asserts that the edges form a forest, so a cycle, multi-edge or self-loop aborts
 *   - Distances are long long sums of edge weights, keep 2 n max|w| below 9.2e18 so d1 + d2 cannot overflow
 *   - Iterative, so deep trees cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct CentroidDecomposition{
    using Branches = vector<vector<pair<int, long long>>>;

    int n;
    vector<tuple<int, int, long long>> edges;
    vector<int> parent, depth;
    /// Nodes renumbered in BFS order with flat adjacency, measured up to 4x faster than per-node vectors at n = 5e5
    vector<int> original, start;
    vector<pair<int, long long>> adj;
    vector<char> removed;
    vector<int> from, size, order;

    CentroidDecomposition(int n) : n(n), parent(n, -1), depth(n, 0) {}

    void add_edge(int u, int v, long long w = 1){
        edges.push_back(make_tuple(u, v, w));
    }

    /// Skips collecting branches, which is about a third of the build time
    void build(){
        decompose([](int) {});
    }

    template <typename Visit>
    void build(Visit visit){
        decompose([&](int c){ visit(original[c], collect_branches(c)); });
    }

    template <typename OnCentroid>
    void decompose(OnCentroid on_centroid){
        renumber();
        removed.assign(n, 0), from.assign(n, -1), size.assign(n, 0);

        /// (any node of a component, centroid of the component it was cut from)
        vector<pair<int, int>> stack;
        for (int s = 0; s < n; s++){
            if (!removed[s]) stack.push_back({s, -1});
            while (!stack.empty()){
                int r, up; tie(r, up) = stack.back();
                stack.pop_back();
                int c = find_centroid(r);
                removed[c] = 1;
                parent[original[c]] = up == -1 ? -1 : original[up];
                depth[original[c]] = up == -1 ? 0 : depth[original[up]] + 1;

                on_centroid(c);
                for (int j = start[c]; j < start[c + 1]; j++){
                    if (!removed[adj[j].first]) stack.push_back({adj[j].first, c});
                }
            }
        }
    }

    void renumber(){
        vector<int> first(n + 1, 0);
        for (auto e : edges) first[get<0>(e) + 1]++, first[get<1>(e) + 1]++;
        for (int u = 0; u < n; u++) first[u + 1] += first[u];
        vector<int> fill(first.begin(), first.end() - 1);
        vector<int> incident(first[n]);
        for (int e = 0; e < (int)edges.size(); e++){
            int u = get<0>(edges[e]), v = get<1>(edges[e]);
            incident[fill[u]++] = e, incident[fill[v]++] = e;
        }

        vector<int> id(n, -1);
        int trees = 0;
        original.clear();
        for (int s = 0; s < n; s++){
            if (id[s] != -1) continue;
            trees++, id[s] = original.size(), original.push_back(s);
            for (int i = id[s]; i < (int)original.size(); i++){
                int u = original[i];
                for (int j = first[u]; j < first[u + 1]; j++){
                    int a = get<0>(edges[incident[j]]), b = get<1>(edges[incident[j]]);
                    int v = a ^ b ^ u;
                    if (id[v] == -1) id[v] = original.size(), original.push_back(v);
                }
            }
        }
        /// A forest is exactly a graph with n - trees edges, so this catches every cycle, multi-edge and self-loop
        assert((int)edges.size() == n - trees);

        start.assign(n + 1, 0), adj.resize(first[n]);
        for (int i = 0; i < n; i++){
            int u = original[i], k = start[i];
            for (int j = first[u]; j < first[u + 1]; j++){
                int a = get<0>(edges[incident[j]]), b = get<1>(edges[incident[j]]); long long w = get<2>(edges[incident[j]]);
                adj[k++] = {id[a ^ b ^ u], w};
            }
            start[i + 1] = k;
        }
    }

    Branches collect_branches(int c){
        Branches branches;

        for (int j = start[c]; j < start[c + 1]; j++){
            int v = adj[j].first; long long w = adj[j].second;
            if (removed[v]) continue;
            vector<pair<int, long long>> branch = {{v, w}};
            from[v] = c;
            for (int i = 0; i < (int)branch.size(); i++){
                int x = branch[i].first; long long d = branch[i].second;  /// a copy: pushing below may reallocate the branch
                for (int k = start[x]; k < start[x + 1]; k++){
                    int y = adj[k].first; long long wy = adj[k].second;
                    if (y != from[x] && !removed[y]) from[y] = x, branch.push_back({y, d + wy});
                }
            }
            for (auto& e : branch) e.first = original[e.first];
            branches.push_back(move(branch));
        }

        return branches;
    }

    int find_centroid(int r){
        order.assign(1, r);
        from[r] = -1, size[r] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j].first;
                if (v == from[u] || removed[v]) continue;
                from[v] = u, size[v] = 0, order.push_back(v);
            }
        }

        /// Nodes holding over half the component form a path down from r, the centroid is its last node in BFS order
        int total = order.size(), centroid = -1;
        for (int i = total - 1; i >= 0; i--){
            int u = order[i];
            size[u]++;
            if (centroid == -1 && 2 * size[u] > total) centroid = u;
            if (i > 0) size[from[u]] += size[u];
        }

        return centroid;
    }
};

long long count_paths_at_most(CentroidDecomposition& cd, long long k){
    auto count_pairs_at_most = [k](vector<long long> d){
        sort(d.begin(), d.end());

        long long pairs = 0;
        int j = (int)d.size() - 1;
        for (int i = 0; i < j; i++){
            while (i < j && d[i] + d[j] > k) j--;
            pairs += j - i;
        }

        return pairs;
    };

    long long pairs = 0;

    cd.build([&](int, const CentroidDecomposition::Branches& branches){
        vector<long long> all = {0};
        for (auto& branch : branches){
            vector<long long> d;
            for (auto e : branch) d.push_back(e.second);
            pairs -= count_pairs_at_most(d);
            all.insert(all.end(), d.begin(), d.end());
        }
        pairs += count_pairs_at_most(all);
    });
    return pairs;
}

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
        for (; isdigit(c); c = read_char()){
            v = v * 10 - (c - '0');
        }
        x = neg ? v : -v;

        return true;
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return 0;

        for (; !isspace(c) && c != EOF; c = read_char()){
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

const int SHIFT = 25;
const long long MASK = (1LL << SHIFT) - 1;

int n, m, res;

void solve(int, const CentroidDecomposition::Branches &branches){
    int i, j, d, l;
    vector<tuple<int, int, int>> items;

    items.push_back(make_tuple(0, 0, -1));
    for (i = 0; i < (int)branches.size(); i++){
        for (j = 0; j < (int)branches[i].size(); j++){
            d = branches[i][j].second >> SHIFT, l = branches[i][j].second & MASK;
            if (d <= m) items.push_back(make_tuple(d, l, i));
        }
    }
    sort(items.begin(), items.end());

    int k = items.size();
    vector<int> dist(k), best1(k), from1(k), best2(k);
    for (i = 0; i < k; i++){
        tie(d, l, j) = items[i];
        dist[i] = d;
        if (i == 0){
            best1[i] = l, from1[i] = j, best2[i] = -1;
            continue;
        }

        best1[i] = best1[i - 1], from1[i] = from1[i - 1], best2[i] = best2[i - 1];
        if (l > best1[i]){
            if (j != from1[i]) best2[i] = best1[i];
            best1[i] = l, from1[i] = j;
        }
        else if (j != from1[i]) best2[i] = max(best2[i], l);
    }

    for (i = 0; i < k; i++){
        tie(d, l, j) = items[i];
        int p = upper_bound(dist.begin(), dist.end(), m - d) - dist.begin() - 1;
        if (p < 0) continue;
        if (from1[p] != j) res = max(res, l + best1[p]);
        else if (best2[p] >= 0) res = max(res, l + best2[p]);
    }
}

int main(){
    int t, i, u, v, d, l, T = 0;
    char buf[64];

    while (read(t)){
        while (t--){
            read(n, m);
            CentroidDecomposition cd(n);
            for (i = 1; i < n; i++){
                read(u, v, d, l);
                cd.add_edge(u - 1, v - 1, ((long long)d << SHIFT) | l);
            }

            res = 0;
            cd.build(solve);

            snprintf(buf, sizeof(buf), "Case %d: %d", ++T, res);
            write(string(buf));
        }
    }

    return 0;
}
