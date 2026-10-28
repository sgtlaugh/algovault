/***
 *
 * Heavy Light Decomposition
 * Maps every tree path to O(log n) contiguous ranges of positions, so any range structure answers path queries
 *
 * Complexity: O(n) to build, O(log n) ranges per path
 *
 * HLD hld(n): an empty tree on nodes 0..n-1
 * hld.add_edge(u, v): adds the edge u - v
 * hld.build(root): roots the tree at root, 0 by default, asserts that every node is reached exactly once
 * hld.update_path(seg, u, v, x, edges = false): seg.update(l, r, x) on every range of the path u - v
 * hld.edge_pos(u, v): the position holding the value of tree edge (u, v), in either order, asserts u and v are adjacent
 * hld.lca(u, v): the lowest common ancestor of u and v
 * hld.ordered_path(u, v, edges = false): the ranges of the path as {from, to} pairs, in path order from u to v
 * hld.path(u, v, edges = false): ranges [l, r] covering the path u - v, in no particular order
 * hld.query_path(seg, u, v, edges = false): the fold of the path through seg.combine, T() for an empty path
 * hld.subtree(v): the single range covering the subtree of v
 * hld.pos[v]: the position of v in [0, n - 1], where the value of node v goes
 * hld.parent[v]: the parent of v, -1 for the root
 * hld.depth[v]: the number of edges from the root to v
 *
 * Notes:
 *   - seg is any size n tree with 0-indexed query(l, r), update(l, r, x) and combine(a, b), like segment_tree.cpp
 *   - edges = true excludes the lca, for values on edges: store the weight of edge (parent[v], v) at pos[v]
 *   - ordered_path walks each range from position from to to: from >= to on the u side, from <= to on the v side
 *   - path and query_path suit commutative folds only (sum, min, max, xor), ordered_path suits any, like affine maps
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct HLD{
    struct Chain{ int head_pos, head_parent; };

    int n;
    vector<pair<int, int>> edge_list;
    vector<int> parent, depth, heavy, head, pos, size;
    vector<Chain> chain;  /// v's head pos and parent in one cache line: lca 2-3x faster on random trees

    HLD(int n) : n(n), parent(n, -1), depth(n), heavy(n, -1), head(n), pos(n), size(n, 1), chain(n) {}

    void add_edge(int u, int v){
        edge_list.push_back({u, v});
    }

    void build(int root = 0){
        size.assign(n, 1), heavy.assign(n, -1), parent.assign(n, -1);  /// build may run again, e.g. with a new root

        /// One flat array instead of a vector per node: the neighbors of u are adj[start[u]..start[u + 1] - 1]
        vector<int> start(n + 1, 0), adj(2 * edge_list.size());
        for (auto [u, v] : edge_list) start[u]++, start[v]++;
        for (int i = 1; i <= n; i++) start[i] += start[i - 1];
        for (auto [u, v] : edge_list) adj[--start[u]] = v, adj[--start[v]] = u;

        vector<int> order = {root};
        parent[root] = -1, depth[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j];
                if (v == parent[u]) continue;
                assert(v != root && parent[v] == -1);  /// a revisit means a cycle, which would grow order forever
                parent[v] = u, depth[v] = depth[u] + 1;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i > 0; i--){
            int v = order[i], p = parent[v];
            size[p] += size[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v;
        }

        /// Heavy child is pushed last so it is popped right after its parent, keeping chains and subtrees contiguous
        vector<int> stack = {root};
        head[root] = root;
        for (int cur = 0; !stack.empty(); cur++){
            int u = stack.back();
            stack.pop_back();
            pos[u] = cur;
            chain[u] = {pos[head[u]], parent[head[u]]};  /// head[u] is u or was popped before u
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j];
                if (v == parent[u] || v == heavy[u]) continue;
                head[v] = v;
                stack.push_back(v);
            }
            if (heavy[u] != -1){
                head[heavy[u]] = head[u];
                stack.push_back(heavy[u]);
            }
        }
    }

    template <typename Tree, typename Value>
    void update_path(Tree& seg, int u, int v, Value x, bool edges = false) const{
        for (auto [l, r] : path(u, v, edges)) seg.update(l, r, x);
    }

    int edge_pos(int u, int v) const{
        assert(parent[u] == v || parent[v] == u);
        return pos[depth[u] > depth[v] ? u : v];  /// edge values live on the child
    }

    /// The head with the larger position is never an ancestor of the other node, so its chain is safe to leave
    int lca(int u, int v) const{
        while (chain[u].head_pos != chain[v].head_pos){
            if (chain[u].head_pos < chain[v].head_pos) swap(u, v);
            u = chain[u].head_parent;
        }
        return pos[u] < pos[v] ? u : v;
    }

    vector<pair<int, int>> ordered_path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> up, down;

        while (chain[u].head_pos != chain[v].head_pos){
            if (chain[u].head_pos > chain[v].head_pos){
                up.push_back({pos[u], chain[u].head_pos});
                u = chain[u].head_parent;
            }
            else{
                down.push_back({chain[v].head_pos, pos[v]});
                v = chain[v].head_parent;
            }
        }

        /// u and v now share a chain, the shallower one (smaller position) is the lca
        if (pos[u] >= pos[v]){
            if (pos[u] >= pos[v] + edges) up.push_back({pos[u], pos[v] + edges});
        }
        else down.push_back({pos[u] + edges, pos[v]});

        up.insert(up.end(), down.rbegin(), down.rend());
        return up;
    }

    vector<pair<int, int>> path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> ranges = ordered_path(u, v, edges);
        for (auto& [a, b] : ranges){
            if (a > b) swap(a, b);
        }
        return ranges;
    }

    /// Folds from the first range, so seg needs no identity; only u == v with edges gives an empty path
    template <typename Tree>
    auto query_path(Tree& seg, int u, int v, bool edges = false) const{
        using T = decltype(seg.query(0, 0));
        auto ranges = path(u, v, edges);
        if (ranges.empty()) return T();

        T res = seg.query(ranges[0].first, ranges[0].second);
        for (int i = 1; i < (int)ranges.size(); i++){
            res = seg.combine(res, seg.query(ranges[i].first, ranges[i].second));
        }
        return res;
    }

    pair<int, int> subtree(int v) const{
        return {pos[v], pos[v] + size[v] - 1};
    }
};

/***
 *
 * Segment Tree with Lazy Propagation
 * Range update, range query and binary search for any associative combine and lazy update
 *
 * Complexity: O(n) to build, O(log n) per update, query, max_right and min_left
 *
 * SegmentTree<T, Ops = AddSum<T>, L = T> st(n): n copies of T(), range add and range sum by default
 * SegmentTree<T, Ops, L> st(ar): built from ar
 * st.update(l, r, lz): apply lz to a[l..r], 0 <= l <= r < n
 * st.max_right(l, pred): largest r with pred(query(l, r)) true, l - 1 if none, 0 <= l <= n
 * st.min_left(r, pred): smallest l with pred(query(l, r)) true, r + 1 if none, -1 <= r < n
 * st.query(l, r): combine over a[l..r], 0 <= l <= r < n
 *
 * Notes:
 *   - pred must be monotone, once false it stays false as the range grows, and is never called on an empty range
 *   - Ops holds three static functions, L is the lazy type when it differs from T, and no identity is needed:
 *       - combine(a, b): fold of two adjacent ranges, a on the left
 *       - apply(val, lz, len): value of a range of len elements after lz
 *       - compose(old, new): one lazy meaning old then new, a lazy that changes nothing must return old
 *   - Arithmetic progression add (x, x + d, ... on [l, r]) with range sum: T = {sum, cnt, idx_sum}, L = {c0, c1}
 *       - leaves are {ar[i], 1, i}, combine and compose add componentwise, apply adds c0 * cnt + c1 * idx_sum
 *       - the update is update(l, r, {x - d * l, d})
 *       - sums reach d * n^2 / 2, so large inputs need a modulus or __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct AddSum{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
};

/// Inherits Ops so seg.combine(a, b) stays callable from outside, hld.cpp's query_path relies on it
template<typename T, typename Ops = AddSum<T>, typename L = T>
struct SegmentTree : Ops{
    int n;
    vector<T> tree;
    vector<L> lazy;
    vector<char> has_lazy;

    SegmentTree(int n) : SegmentTree(vector<T>(n)) {}

    SegmentTree(const vector<T>& ar) : n(ar.size()), tree(4 * n), lazy(4 * n), has_lazy(4 * n){
        if (n) build(ar, 1, 0, n - 1);
    }

    void update(int l, int r, L val){
        assert(0 <= l && l <= r && r < n);
        update(1, 0, n - 1, l, r, val);
    }

    void update(int idx, int a, int b, int l, int r, const L& val){
        if (a == l && b == r){
            push_lazy(idx, val);
            propagate(idx, a, b);
            return;
        }
        propagate(idx, a, b);

        /// the untouched child gets propagated too, so the combine below reads its true value
        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c){
            propagate(q, c + 1, b);
            update(p, a, c, l, r, val);
        }
        else if (l > c){
            propagate(p, a, c);
            update(q, c + 1, b, l, r, val);
        }
        else{
            update(p, a, c, l, c, val);
            update(q, c + 1, b, c + 1, r, val);
        }
        tree[idx] = Ops::combine(tree[p], tree[q]);
    }

    template<typename Pred>
    int max_right(int l, Pred pred){
        assert(0 <= l && l <= n);
        T acc{};
        bool empty = true;
        int pos = l < n ? first_fail(1, 0, n - 1, l, pred, acc, empty) : -1;
        return pos == -1 ? n - 1 : pos - 1;
    }

    template<typename Pred>
    int min_left(int r, Pred pred){
        assert(-1 <= r && r < n);
        T acc{};
        bool empty = true;
        int pos = r >= 0 ? last_fail(1, 0, n - 1, r, pred, acc, empty) : -1;
        return pos == -1 ? 0 : pos + 1;
    }

    T query(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r);
    }

    T query(int idx, int a, int b, int l, int r){
        propagate(idx, a, b);
        if (a == l && b == r) return tree[idx];

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c) return query(p, a, c, l, r);
        if (l > c) return query(q, c + 1, b, l, r);
        return Ops::combine(query(p, a, c, l, c), query(q, c + 1, b, c + 1, r));
    }

    void build(const vector<T>& ar, int idx, int a, int b){
        if (a == b){
            tree[idx] = ar[a];
            return;
        }

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        build(ar, p, a, c);
        build(ar, q, c + 1, b);
        tree[idx] = Ops::combine(tree[p], tree[q]);
    }

    void push_lazy(int idx, const L& val){
        lazy[idx] = has_lazy[idx] ? Ops::compose(lazy[idx], val) : val;
        has_lazy[idx] = 1;
    }

    void propagate(int idx, int a, int b){
        if (!has_lazy[idx]) return;
        tree[idx] = Ops::apply(tree[idx], lazy[idx], b - a + 1);
        if (a != b){
            push_lazy(idx << 1, lazy[idx]);
            push_lazy(idx << 1 | 1, lazy[idx]);
        }
        has_lazy[idx] = 0;
    }

    /// First position in [max(a, l), b] where pred(acc merged with l..pos) fails, or -1, acc absorbs each passing node
    template<typename Pred>
    int first_fail(int idx, int a, int b, int l, Pred& pred, T& acc, bool& empty){
        if (b < l) return -1;
        propagate(idx, a, b);
        if (l <= a){
            T nxt = empty ? tree[idx] : Ops::combine(acc, tree[idx]);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = first_fail(idx << 1, a, c, l, pred, acc, empty);
        return res != -1 ? res : first_fail(idx << 1 | 1, c + 1, b, l, pred, acc, empty);
    }

    /// mirror of first_fail scanning leftwards from r, combine order kept as values pos..r for non-commutative combines
    template<typename Pred>
    int last_fail(int idx, int a, int b, int r, Pred& pred, T& acc, bool& empty){
        if (a > r) return -1;
        propagate(idx, a, b);
        if (b <= r){
            T nxt = empty ? tree[idx] : Ops::combine(tree[idx], acc);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = last_fail(idx << 1 | 1, c + 1, b, r, pred, acc, empty);
        return res != -1 ? res : last_fail(idx << 1, a, c, r, pred, acc, empty);
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

struct Ops{
    static int combine(int a, int b){ return a + b; }
    static int apply(int val, int lz, int len){ return lz ? len - val : val; }
    static int compose(int old_lz, int new_lz){ return old_lz ^ new_lz; }
};

int n, k, deg[MAX], idx[MAX], root_of[MAX], eu[MAX], ev[MAX];
vector<int> adj[MAX], cyc;

void find_cycle(){
    int i, u;
    queue<int> Q;

    for (i = 0; i < n; i++){
        if (deg[i] == 1) Q.push(i);
    }
    while (!Q.empty()){
        u = Q.front();
        Q.pop();
        for (int v : adj[u]){
            if (--deg[v] == 1) Q.push(v);
        }
    }

    for (i = 0; deg[i] < 2; i++){}
    int prev = -1, cur = i;
    do{
        idx[cur] = cyc.size();
        cyc.push_back(cur);
        for (int v : adj[cur]){
            if (deg[v] >= 2 && v != prev && v != cyc[0]){
                prev = cur, cur = v;
                break;
            }
        }
    } while (cur != cyc.back());
    k = cyc.size();
}

void find_roots(){
    vector<int> st = cyc;

    memset(root_of, -1, sizeof(int) * n);
    for (int c : cyc) root_of[c] = c;
    while (!st.empty()){
        int u = st.back();
        st.pop_back();
        for (int v : adj[u]){
            if (root_of[v] == -1){
                root_of[v] = root_of[u];
                st.push_back(v);
            }
        }
    }
}

int main(){
    int m, i, a, b, extra = 0;

    read(n, m);
    for (i = 0; i < n; i++){
        read(eu[i], ev[i]);
        eu[i]--, ev[i]--;
        adj[eu[i]].push_back(ev[i]);
        adj[ev[i]].push_back(eu[i]);
        deg[eu[i]]++, deg[ev[i]]++;
    }

    find_cycle();
    find_roots();

    HLD hld(n);
    for (i = 0; i < n; i++){
        if (min(eu[i], ev[i]) == min(cyc[0], cyc[k - 1]) && max(eu[i], ev[i]) == max(cyc[0], cyc[k - 1])) continue;
        hld.add_edge(eu[i], ev[i]);
    }
    hld.build(cyc[0]);
    SegmentTree<int, Ops> seg(n);

    while (m--){
        read(a, b);
        a--, b--;

        int x = idx[root_of[a]], y = idx[root_of[b]];
        int d = abs(x - y), dir = (y > x) ? 1 : -1;
        if (root_of[a] == root_of[b] || (2 * d) < k || ((2 * d) == k && cyc[x + dir] < cyc[(x - dir + k) % k])){
            hld.update_path(seg, a, b, 1, true);
        }
        else{
            hld.update_path(seg, a, cyc[x < y ? 0 : k - 1], 1, true);
            hld.update_path(seg, cyc[x < y ? k - 1 : 0], b, 1, true);
            extra ^= 1;
        }

        int on = seg.query(0, n - 1) + extra;
        int cyc_on = hld.query_path(seg, cyc[0], cyc[k - 1], true) + extra;
        write(n - on + (cyc_on == k));
    }
    return 0;
}
