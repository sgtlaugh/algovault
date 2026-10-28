/***
 *
 * Rerooting DP
 * A tree DP evaluated with every node as the root, from three user operations
 *
 * Complexity: O(n), 2(n - 1) calls to apply_edge and O(n) calls to combine
 *
 * Rerooting<T, E = long long> tree(n): nodes 0..n-1, T is the DP value of a rooted subtree and E the edge data
 * tree.add_edge(u, v, w = E()): adds the edge (u, v) with data w, 0 <= u, v < n
 * tree.solve(combine, identity, apply_edge): res with res[r] the contributions of every neighbor of r combined
 *
 * Notes:
 *   - apply_edge(value, from, to, w): what the subtree at from, whose DP is value, contributes to to across edge w
 *   - combine(a, b): two contributions as one, associative and commutative with combine(identity, a) = a
 *   - solve asserts n >= 1 and a tree, res[0] is identity when n = 1
 *   - Node data goes inside apply_edge through from, the child's own node, or is applied to res afterwards
 *   - E must be default constructible
 *   - Iterative, so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename E = long long>
struct Rerooting{
    int n;
    vector<array<int, 2>> ends;
    vector<E> data;

    Rerooting(int n) : n(n) {}

    void add_edge(int u, int v, E w = E()){
        assert(0 <= u && u < n && 0 <= v && v < n);
        ends.push_back({u, v});
        data.push_back(w);
    }

    template<typename Combine, typename Apply>
    vector<T> solve(Combine combine, const T& identity, Apply apply_edge) const{
        assert(n >= 1 && (int)ends.size() == n - 1);

        /// The edges at v are adj[start[v]..start[v + 1] - 1] in insertion order, as indices into ends and data
        vector<int> start(n + 1, 0), adj(2 * (n - 1));
        for (auto [u, v] : ends) start[u + 1]++, start[v + 1]++;
        for (int v = 0; v < n; v++) start[v + 1] += start[v];
        vector<int> pos(start.begin(), start.end() - 1);
        for (int e = 0; e < n - 1; e++) adj[pos[ends[e][0]]++] = e, adj[pos[ends[e][1]]++] = e;

        /// Indexed by BFS position from 0, so the children of i are first[i]..first[i + 1] - 1 and passes scan in order
        vector<int> order = {0}, parent(n, -1), edge = {-1}, first(n + 1, n);
        vector<char> seen(n, false);
        seen[0] = true;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            first[i] = order.size();
            for (int j = start[u]; j < start[u + 1]; j++){
                int e = adj[j], to = ends[e][0] ^ ends[e][1] ^ u;
                if (e == edge[i]) continue;
                assert(!seen[to]);
                seen[to] = true, parent[order.size()] = i;
                order.push_back(to), edge.push_back(e);
            }
        }
        assert((int)order.size() == n);

        /// down[i] combines i's children, then what i's subtree gives its parent, up[i] what the parent side gives i
        vector<T> down(n, identity), up(n, identity), suffix, res(n, identity);
        for (int i = n - 1; i > 0; i--){
            down[i] = apply_edge(down[i], order[i], order[parent[i]], data[edge[i]]);
            down[parent[i]] = combine(down[parent[i]], down[i]);
        }

        for (int i = 0; i < n; i++){
            int l = first[i], k = first[i + 1] - l;
            suffix.assign(k + 1, identity);
            for (int j = k - 1; j >= 0; j--) suffix[j] = combine(down[l + j], suffix[j + 1]);
            res[order[i]] = combine(up[i], suffix[0]);

            T prefix = up[i];
            for (int j = 0; j < k; j++){
                up[l + j] = apply_edge(combine(prefix, suffix[j + 1]), order[i], order[l + j], data[edge[l + j]]);
                prefix = combine(prefix, down[l + j]);
            }
        }

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

#define MAX 100010

int n, len[MAX], parent[MAX];
bool is_file[MAX];

pair<long long, long long> lift(pair<long long, long long> x, int from, int to, long long){
    long long cnt = x.first + is_file[from];
    if (parent[from] == to) return make_pair(cnt, x.second + cnt * (len[from] + 1));
    return make_pair(cnt, x.second + cnt * 3);
}

int main(){
    int i, j, m, v;
    string name;

    freopen("dirtraverse.in", "r", stdin);
    freopen("dirtraverse.out", "w", stdout);

    read(n);
    Rerooting<pair<long long, long long>> tree(n);
    for (i = 0; i < n; i++){
        read(name, m);
        len[i] = name.size(), is_file[i] = (m == 0);
        for (j = 0; j < m; j++){
            read(v);
            parent[v - 1] = i;
            tree.add_edge(i, v - 1);
        }
    }

    auto add = [](pair<long long, long long> a, pair<long long, long long> b){
        return make_pair(a.first + b.first, a.second + b.second);
    };
    vector<pair<long long, long long>> sums = tree.solve(add, make_pair(0LL, 0LL), lift);

    long long res = LLONG_MAX;
    for (i = 0; i < n; i++){
        if (!is_file[i]) res = min(res, sums[i].second - sums[i].first);
    }

    write(res);
    return 0;
}
