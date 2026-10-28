/***
 *
 * Bridges and the Bridge Tree
 * Finds the bridges and 2-edge-connected components of an undirected graph, parallel edges and self loops allowed
 *
 * Complexity: O(n + m)
 *
 * BridgeGraph g(n): nodes are 0 to n - 1
 * g.add_edge(u, v): adds an edge and returns its id, ids count up from 0 in add_edge order
 * g.get_bridge_tree(): adjacency lists over the g.count components, one tree edge per bridge, a forest
 * g.get_bridges(): each edge whose removal adds a component, with endpoints u, v, id and node counts cnt_u, cnt_v
 * g.get_components(): the 2-edge-connected components, the groups left after removing every bridge, each ascending
 * g.comp[v]: component of v in [0, g.count), numbered by smallest vertex, set by get_components and get_bridge_tree
 *
 * Notes:
 *   - A doubled edge is never a bridge
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Bridge{
    int u, v;
    int cnt_u, cnt_v; /// node counts on the side of u and of v once the bridge is removed
    int id;
};

struct BridgeGraph{
    int n, m = 0, count = 0;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge index)
    vector<int> comp;

    BridgeGraph(int n): n(n), adj(n), comp(n, -1) {}

    int add_edge(int u, int v){
        adj[u].push_back({v, m});
        adj[v].push_back({u, m});
        return m++;
    }

    vector<vector<int>> get_bridge_tree(){
        get_components();

        vector<vector<int>> tree(count);
        for (int u = 0; u < n; u++){
            for (auto e : adj[u]){
                int v = e.first;
                if (comp[u] != comp[v]) tree[comp[u]].push_back(comp[v]);
            }
        }

        return tree;
    }

    /// Bridges are listed in the order the DFS finishes their endpoint v
    vector<Bridge> get_bridges(){
        vector<Bridge> bridges;
        vector<int> discover(n, 0), low(n), parent_edge(n), edge_pos(n, 0), dfs_stack;

        for (int i = 0; i < n; i++){
            if (discover[i]) continue;

            int first = bridges.size(), dt = 0;
            discover[i] = low[i] = ++dt, parent_edge[i] = -1;
            dfs_stack.push_back(i);

            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                if (edge_pos[u] < (int)adj[u].size()){
                    int v, id;
                    tie(v, id) = adj[u][edge_pos[u]++];
                    if (!discover[v]){
                        discover[v] = low[v] = ++dt, parent_edge[v] = id;
                        dfs_stack.push_back(v);
                    }
                    else if (id != parent_edge[u]) low[u] = min(low[u], discover[v]);
                    continue;
                }

                dfs_stack.pop_back();
                if (dfs_stack.empty()) break;

                int p = dfs_stack.back();
                low[p] = min(low[p], low[u]);
                if (low[u] > discover[p]) bridges.push_back({p, u, 0, dt - discover[u] + 1, parent_edge[u]});
            }

            /// The component size is known only once its DFS ends
            for (int j = first; j < (int)bridges.size(); j++) bridges[j].cnt_u = dt - bridges[j].cnt_v;
        }

        return bridges;
    }

    vector<vector<int>> get_components(){
        vector<char> is_bridge(m, 0);
        for (auto& b : get_bridges()) is_bridge[b.id] = 1;

        count = 0;
        vector<int> dfs_stack;
        fill(comp.begin(), comp.end(), -1);
        for (int root = 0; root < n; root++){
            if (comp[root] != -1) continue;
            comp[root] = count;
            dfs_stack.push_back(root);
            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                dfs_stack.pop_back();
                for (auto e : adj[u]){
                    int v = e.first, id = e.second;
                    if (comp[v] == -1 && !is_bridge[id]) comp[v] = count, dfs_stack.push_back(v);
                }
            }
            count++;
        }

        vector<vector<int>> groups(count);
        for (int v = 0; v < n; v++) groups[comp[v]].push_back(v);
        return groups;
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

int n, m, A[MAX], B[MAX];
string name[20010];
unordered_map<string, int> mp;

int get_id(const string& s){
    if (!mp.count(s)){
        int id = mp.size();
        mp[s] = id;
        name[id] = s;
    }

    return mp[s];
}

int main(){
    int t, i;
    string a, b;

    while (read(t)){
        while (t--){
            read(n, m);
            mp.clear();
            for (i = 0; i < m; i++){
                read(a, b);
                A[i] = get_id(a), B[i] = get_id(b);
            }

            if ((int)mp.size() < n){
                write(0);
                continue;
            }

            BridgeGraph g(n);
            for (i = 0; i < m; i++) g.add_edge(A[i], B[i]);
            g.get_components();
            vector<Bridge> bridges = g.get_bridges();
            if ((int)bridges.size() != (g.count - 1)){
                write(0);
                continue;
            }

            vector<int> ids;
            for (i = 0; i < (int)bridges.size(); i++) ids.push_back(bridges[i].id);
            sort(ids.begin(), ids.end());

            write((int)ids.size());
            for (i = 0; i < (int)ids.size(); i++) write(name[A[ids[i]]], name[B[ids[i]]]);
        }
    }

    return 0;
}
