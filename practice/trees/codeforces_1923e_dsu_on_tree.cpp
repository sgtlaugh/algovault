/***
 *
 * DSU on Tree (small to large, sack)
 * Answers a query about every subtree by keeping the heavy child's data and re-adding only light subtrees
 *
 * Complexity: O(n log n) calls to add / remove
 *
 * DsuOnTree tree(n): an empty tree on nodes 0..n-1
 * tree.add_edge(u, v): adds the edge u - v
 * tree.run(root, add, remove, answer): calls answer(v) for every v while the structure holds exactly v's subtree
 * tree.tin[v]: the position of v in order
 * tree.tout[v]: the last position of v's subtree, so the subtree of v is order[tin[v] .. tout[v]]
 * tree.order[i]: the node at position i of a preorder that visits each node's heavy child first
 * tree.parent[v]: the parent of v, -1 for the root
 * tree.size[v]: the number of nodes in the subtree of v
 * tree.heavy[v]: a child of v with the largest subtree, -1 for a leaf
 *
 * Notes:
 *   - The fields describe the tree rooted at the root of the last run call
 *   - add(v) / remove(v) put node v into / take it out of your data structure, empty again when run returns
 *   - The graph must be a tree, run asserts that every node is reached exactly once
 *   - Iterative, so deep trees cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DsuOnTree{
    int n;
    vector<vector<int>> adj;
    vector<int> parent, heavy, size, tin, tout, order;

    DsuOnTree(int n) : n(n), adj(n), parent(n, -1), heavy(n, -1), size(n, 1), tin(n), tout(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    template <typename Add, typename Remove, typename Answer>
    void run(int root, Add add, Remove remove, Answer answer){
        prepare(root);

        /// Reverse preorder reaches v right after its heavy child, whose subtree the structure still holds
        for (int i = n - 1; i >= 0; i--){
            int v = order[i];
            add(v);
            for (int j = heavy[v] == -1 ? i + 1 : tout[heavy[v]] + 1; j <= tout[v]; j++) add(order[j]);
            answer(v);
            if (v == root || heavy[parent[v]] != v){
                for (int j = i; j <= tout[v]; j++) remove(order[j]);
            }
        }
    }

    /// Preorder with the heavy child first, so every subtree is one contiguous range
    void prepare(int root){
        size.assign(n, 1), heavy.assign(n, -1), parent.assign(n, -1);  /// run may be called again, from any root
        vector<int> bfs = {root};
        parent[root] = -1;
        for (int i = 0; i < (int)bfs.size(); i++){
            for (int v : adj[bfs[i]]){
                if (v == parent[bfs[i]]) continue;
                assert(v != root && parent[v] == -1);  /// a revisit means a cycle, which would grow bfs forever
                parent[v] = bfs[i], bfs.push_back(v);
            }
        }
        assert((int)bfs.size() == n);

        for (int i = n - 1; i > 0; i--){
            int v = bfs[i], p = parent[v];
            size[p] += size[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v;
        }

        order.clear();
        vector<int> stack = {root};
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            tin[u] = order.size(), order.push_back(u);
            for (int v : adj[u]){
                if (v != parent[u] && v != heavy[u]) stack.push_back(v);
            }
            if (heavy[u] != -1) stack.push_back(heavy[u]);
        }

        for (int i = n - 1; i >= 0; i--) tout[order[i]] = tin[order[i]] + size[order[i]] - 1;
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

#define MAX 200010

int n, color[MAX], cnt[MAX], path_cnt[MAX];
bool present[MAX];
long long res;
vector<int> path;
DsuOnTree *tree;

void clear_path(){
    while (!path.empty()){
        path_cnt[color[path.back()]]--;
        path.pop_back();
    }
}

void add_node(int v){
    int c = color[v], p = tree->parent[v];

    if (p == -1 || !present[p]){
        clear_path();
        cnt[c] = 1;
    }
    else{
        while (path.back() != p){
            path_cnt[color[path.back()]]--;
            path.pop_back();
        }
        if (path_cnt[c] == 0) cnt[c]++;
    }

    present[v] = true;
    path.push_back(v);
    path_cnt[c]++;
}

void remove_node(int v){
    present[v] = false;
    cnt[color[v]] = 0;
}

void answer(int v){
    if (tree->parent[v] == -1){
        for (int c = 1; c <= n; c++) res += (long long)cnt[c] * (cnt[c] - 1) / 2;
        return;
    }

    long long s = cnt[color[tree->parent[v]]];
    res += (s * (s - 1) / 2) + s;
}

int main(){
    int t, i, u, v;

    read(t);
    while (t--){
        read(n);
        for (i = 0; i < n; i++) read(color[i]);

        DsuOnTree cur(n);
        for (i = 1; i < n; i++){
            read(u, v);
            cur.add_edge(u - 1, v - 1);
        }

        res = 0, tree = &cur;
        cur.run(0, add_node, remove_node, answer);
        clear_path();
        write(res);
    }
    return 0;
}
