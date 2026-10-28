/***
 *
 * Dominator Tree (Lengauer-Tarjan)
 * Immediate dominator of every vertex of a directed graph, for a fixed root
 *
 * Complexity: O((n + m) log n)
 *
 * DominatorTree g(n): vertices are 0 to n - 1, parallel edges and self loops allowed
 * g.add_edge(u, v): directed edge u -> v
 * g.run(root): idom, idom[v] the closest strict dominator of v, idom[root] = root, -1 if v is unreachable from root
 *
 * Notes:
 *   - d dominates v if every path from root to v passes through d
 *   - idom is the dominator tree as a parent array: d dominates v iff d is an ancestor of v in it
 *   - run can be called again with another root
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DominatorTree{
    int n;
    vector<vector<int>> adj;
    vector<int> order, num, parent, sdom, label, anc, path;

    DominatorTree(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    vector<int> run(int root){
        dfs(root);
        int k = order.size();
        vector<vector<int>> pred(k), bucket(k);
        for (int i = 0; i < k; i++){
            for (int v : adj[order[i]]) pred[num[v]].push_back(i);
        }

        vector<int> dom(k, 0);
        sdom.resize(k), label.resize(k), anc.resize(k);
        for (int i = 0; i < k; i++) sdom[i] = label[i] = anc[i] = i;
        for (int i = k - 1; i >= 0; i--){
            for (int p : pred[i]) sdom[i] = min(sdom[i], sdom[eval(p)]);
            if (i > 0) bucket[sdom[i]].push_back(i);
            for (int w : bucket[i]){
                int v = eval(w);
                dom[w] = sdom[v] == sdom[w] ? sdom[w] : v;
            }
            if (i > 0) anc[i] = parent[i];
        }

        vector<int> idom(n, -1);
        idom[root] = root;
        for (int i = 1; i < k; i++){
            if (dom[i] != sdom[i]) dom[i] = dom[dom[i]];
            idom[order[i]] = order[dom[i]];
        }

        return idom;
    }

    /// Preorder numbering, parent[i] is the DFS tree parent of the i-th visited vertex, num[v] = -1 if unreachable
    void dfs(int root){
        order.clear(), parent.clear();
        num.assign(n, -1);
        vector<int> it(n, 0), stack = {root};
        num[root] = 0, order.push_back(root), parent.push_back(0);
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == (int)adj[u].size()){
                stack.pop_back();
                continue;
            }

            int v = adj[u][it[u]++];
            if (num[v] != -1) continue;
            num[v] = order.size();
            order.push_back(v), parent.push_back(num[u]);
            stack.push_back(v);
        }
    }

    /// Vertex of minimum sdom on the forest path from v up to its root, the root itself excluded
    int eval(int v){
        if (anc[v] == v) return v;

        path.clear();
        for (int u = v; anc[anc[u]] != anc[u]; u = anc[u]) path.push_back(u);
        for (int i = (int)path.size() - 1; i >= 0; i--){
            int w = path[i], a = anc[w];
            if (sdom[label[a]] < sdom[label[w]]) label[w] = label[a];
            anc[w] = anc[a];
        }

        return label[v];
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

#define MAX 100010
#define LOG 17

int n, m, t, depth[MAX], tin[MAX], tout[MAX], up[LOG][MAX], mark[MAX];
unsigned long long reach_from[MAX], reach_to[MAX], below[MAX], bad[MAX];
vector<int> adj[MAX], radj[MAX], child[MAX];

int lca(int u, int v){
    if (depth[u] < depth[v]) swap(u, v);
    for (int j = LOG - 1; j >= 0; j--){
        if ((depth[u] - (1 << j)) >= depth[v]) u = up[j][u];
    }
    if (u == v) return u;
    for (int j = LOG - 1; j >= 0; j--){
        if (up[j][u] != up[j][v]) u = up[j][u], v = up[j][v];
    }

    return up[0][u];
}

int deepest(int c, int b){
    if ((reach_from[c] >> b) & 1) return depth[c];
    for (int j = LOG - 1; j >= 0; j--){
        int x = up[j][c];
        if (!((reach_from[x] >> b) & 1)) c = x;
    }

    return depth[c] - 1;
}

long long free_nodes(int p, int c, unsigned long long r){
    int lo = depth[p], hi = depth[c];

    if ((hi - lo) <= 1) return 0;

    vector<pair<int, int>> seg;
    for (int b = 0; b < 64; b++){
        if (!((r >> b) & 1)) continue;
        int x = max(lo, depth[lca(c, mark[b])]), y = min(hi - 1, deepest(c, b));
        if (x < y) seg.push_back({x, y});
    }
    sort(seg.begin(), seg.end());

    long long res = hi - lo - 1, cur = lo;
    for (auto& s: seg){
        if (s.second <= cur) continue;
        res -= s.second - max((long long)s.first, cur);
        cur = s.second;
    }

    return res;
}

int main(){
    int q, i, j, k, a, x, y, w = 0;

    read(n, m, q);
    DominatorTree g(n);
    vector<int> indeg(n, 0);
    for (i = 0; i < m; i++){
        read(x, y);
        x--, y--;
        g.add_edge(x, y);
        adj[x].push_back(y), radj[y].push_back(x), indeg[y]++;
    }
    vector<int> idom = g.run(0);

    vector<int> order, deg = indeg;
    for (i = 0; i < n; i++){
        if (!deg[i]) order.push_back(i);
    }
    for (i = 0; i < (int)order.size(); i++){
        for (int v: adj[order[i]]){
            if (!--deg[v]) order.push_back(v);
        }
    }

    for (i = 0; i < n; i++){
        if (indeg[i] >= 2) mark[w] = i, reach_from[i] = reach_to[i] = below[i] = 1ULL << w, w++;
    }
    for (int u: order){
        for (int v: radj[u]) reach_to[u] |= reach_to[v];
    }
    for (i = n - 1; i >= 0; i--){
        int u = order[i];
        for (int v: adj[u]) reach_from[u] |= reach_from[v];
    }

    for (i = 1; i < n; i++) child[idom[i]].push_back(i);
    vector<int> st = {0}, pre;
    up[0][0] = 0;
    while (!st.empty()){
        int u = st.back();
        st.pop_back();
        tin[u] = t++, pre.push_back(u);
        for (int v: child[u]) depth[v] = depth[u] + 1, up[0][v] = u, st.push_back(v);
    }
    for (i = n - 1; i >= 0; i--){
        int u = pre[i];
        tout[u] = tin[u];
        for (int v: child[u]) tout[u] = max(tout[u], tout[v]), below[u] |= below[v];
        bad[u] = reach_from[u] & ~below[u];
    }
    for (j = 1; j < LOG; j++){
        for (i = 0; i < n; i++) up[j][i] = up[j - 1][up[j - 1][i]];
    }

    vector<int> in_a(n, -1), vparent(n);
    vector<array<long long, 5>> f(n);
    for (int qi = 0; qi < q; qi++){
        read(k, a);
        vector<int> v(a);
        read(v);

        unsigned long long r = 0;
        for (auto& z: v) z--, in_a[z] = qi, r |= reach_to[z];

        v.push_back(0);
        sort(v.begin(), v.end(), [&](int p1, int p2){ return tin[p1] < tin[p2]; });
        int s = v.size();
        for (i = 0; i + 1 < s; i++) v.push_back(lca(v[i], v[i + 1]));
        sort(v.begin(), v.end(), [&](int p1, int p2){ return tin[p1] < tin[p2]; });
        v.erase(unique(v.begin(), v.end()), v.end());

        vector<int> stk;
        for (int u: v){
            while (!stk.empty() && tout[stk.back()] < tin[u]) stk.pop_back();
            vparent[u] = stk.empty() ? -1 : stk.back();
            stk.push_back(u);
        }

        for (int u: v){
            f[u].fill(0);
            f[u][0] = (in_a[u] == qi) ? 0 : 1;
        }
        for (i = v.size() - 1; i >= 0; i--){
            int u = v[i];
            if (!(bad[u] & r)) f[u][1]++;
            if (vparent[u] == -1) continue;

            int p = vparent[u];
            array<long long, 5> h = f[u], res = {0, 0, 0, 0, 0};
            h[1] += free_nodes(p, u, r);
            for (x = 0; x <= k; x++){
                for (y = 0; x + y <= k; y++) res[x + y] += f[p][x] * h[y];
            }
            f[p] = res;
        }

        write(f[0][k]);
    }

    return 0;
}
