/***
 *
 * K Shortest Walks (Eppstein)
 * Lengths of the k shortest walks from s to t in a directed graph, vertices and edges may repeat
 *
 * Complexity: O((n + m) log m + k log k) time, O((n + m) log m + k) memory
 *
 * KShortestWalks g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v, w): directed edge, w >= 0
 * g.solve(s, t, k): the k smallest walk lengths in non-decreasing order, fewer if fewer walks exist
 *
 * Notes:
 *   - Lengths count with multiplicity and parallel edges give distinct walks
 *   - When s == t the empty walk counts, so the first value is 0
 *   - The largest reported length plus n max(w) must fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY LeftistHeap from code_library/data_structures/leftist_heap.cpp
template <typename T, typename Compare = less<T>>
struct LeftistHeap{
    static constexpr int EMPTY = -1;

    struct Node{
        T val;
        int left, right, rank;
    };

    vector<Node> nodes;
    Compare cmp;

    /// Copies only the nodes on the merged right spines, whose length is bounded by the ranks
    int meld(int a, int b){
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        if (cmp(nodes[b].val, nodes[a].val)) swap(a, b);

        Node copy = nodes[a];
        nodes.push_back(copy);
        int c = (int)nodes.size() - 1;

        int right = meld(copy.right, b);  /// the call can reallocate nodes, assign after it
        nodes[c].right = right;
        if (rank(nodes[c].left) < rank(right)) swap(nodes[c].left, nodes[c].right);
        nodes[c].rank = rank(nodes[c].right) + 1;
        return c;
    }

    int pop(int h){
        assert(h != EMPTY);
        return meld(nodes[h].left, nodes[h].right);
    }

    int push(int h, const T& x){
        nodes.push_back({x, EMPTY, EMPTY, 1});
        return meld(h, (int)nodes.size() - 1);
    }

    T top(int h) const{
        assert(h != EMPTY);
        return nodes[h].val;
    }

    int rank(int h) const{
        return h == EMPTY ? 0 : nodes[h].rank;
    }
};
/// END COPY LeftistHeap

template <typename T, typename Compare> constexpr int LeftistHeap<T, Compare>::EMPTY;

struct KShortestWalks{
    static constexpr long long INF = numeric_limits<long long>::max();

    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;
    LeftistHeap<pair<long long, int>> pool;

    KShortestWalks(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    vector<long long> solve(int s, int t, int k){
        vector<vector<int>> rev(n), out(n);
        for (int id = 0; id < (int)edges.size(); id++){
            rev[edges[id].v].push_back(id);
            out[edges[id].u].push_back(id);
        }

        vector<long long> d(n, INF);
        vector<int> tree_edge(n, -1), order;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> heap;
        d[t] = 0, heap.push({0, t});
        while (!heap.empty()){
            long long du = heap.top().first;
            int u = heap.top().second;
            heap.pop();
            if (du != d[u]) continue;
            order.push_back(u);
            for (int id : rev[u]){
                int x = edges[id].u;
                long long w = edges[id].w;
                if (du + w < d[x]) d[x] = du + w, tree_edge[x] = id, heap.push({d[x], x});
            }
        }

        if (d[s] == INF || k <= 0) return {};

        /// Dijkstra pops a vertex after its tree parent, so h[parent] is complete when h[u] needs it
        pool.nodes.clear();
        vector<int> h(n, pool.EMPTY);
        for (int u : order){
            vector<pair<long long, int>> sidetracks;
            for (int id : out[u]){
                int v = edges[id].v;
                long long w = edges[id].w;
                if (id != tree_edge[u] && d[v] != INF) sidetracks.push_back({w + d[v] - d[u], v});
            }
            sort(sidetracks.rbegin(), sidetracks.rend());

            /// Sorted decreasing, each push of a smaller value becomes the root above the old heap, copying O(1) nodes
            int own = pool.EMPTY;
            for (auto sidetrack : sidetracks) own = pool.push(own, sidetrack);
            h[u] = pool.meld(own, tree_edge[u] == -1 ? pool.EMPTY : h[edges[tree_edge[u]].v]);
        }

        vector<long long> result = {d[s]};
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> best;
        if (h[s] != pool.EMPTY) best.push({d[s] + pool.top(h[s]).first, h[s]});
        while ((int)result.size() < k && !best.empty()){
            long long len = best.top().first;
            int x = best.top().second;
            best.pop();
            result.push_back(len);

            const auto& node = pool.nodes[x];
            long long key = node.val.first;
            int to = node.val.second;
            if (node.left != pool.EMPTY) best.push({len - key + pool.top(node.left).first, node.left});
            if (node.right != pool.EMPTY) best.push({len - key + pool.top(node.right).first, node.right});
            if (h[to] != pool.EMPTY) best.push({len + pool.top(h[to]).first, h[to]});
        }

        return result;
    }
};

constexpr long long KShortestWalks::INF;

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
    int n, m, i, u, v, w, s, t, k;

    while (read(n, m)){
        KShortestWalks g(n);
        for (i = 0; i < m; i++){
            read(u, v, w);
            g.add_edge(u - 1, v - 1, w);
        }
        read(s, t, k);
        if (s == t) k++;

        vector<long long> res = g.solve(s - 1, t - 1, k);
        write((int)res.size() < k ? -1 : res[k - 1]);
    }

    return 0;
}
