/***
 *
 * Weighted General Graph Matching
 * Matching of maximum total edge weight in a general graph, not necessarily perfect or of maximum size
 *
 * Complexity: O(n^3), O(n^2) memory
 *
 * WeightedBlossom wb(n): vertices 0 to n - 1, no edges
 * wb.add_edge(u, v, w): undirected edge of weight w, u != v, a repeated edge keeps its largest weight
 * wb.maximum_matching(): the largest total weight of a matching, callable again after adding edges
 * wb.mate[v]: after maximum_matching(), the vertex matched to v or -1
 *
 * Notes:
 *   - Weights are integers up to 1e18 / n, an edge with w <= 0 never helps and is dropped
 *   - Heaviest among the largest matchings: add B = n max w / 2 + 1 to every weight, B + max w within the bound above
 *   - Reimplemented from Z. Galil, Efficient Algorithms for Finding Maximum Matching in Graphs, ACM Comput. Surv. 1986
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct WeightedBlossom{
    enum Label{UNLABELED, OUTER, INNER};
    static constexpr pair<int, int> NONE = {-1, -1};

    int n;
    vector<vector<long long>> weight;
    vector<int> mate;

    /// Ids 0 to n - 1 are vertices, n to 2n - 1 blossoms; base[b] == -1 marks an unused blossom id
    vector<long long> dual;
    vector<int> top, parent, base, label, unused_ids, queue;
    /// label_edge[b] = (outer vertex outside b, vertex inside b), NONE.first for a tree root
    vector<pair<int, int>> label_edge, best_edge;
    /// child_edges[b][i] joins children[b][i] to the next child, children[b][0] holds the base
    vector<vector<int>> children;
    vector<vector<pair<int, int>>> child_edges, best_edges_to;
    vector<char> built_this_stage, on_path;

    WeightedBlossom(int n): n(n), weight(n, vector<long long>(n, 0)), mate(n, -1) {}

    void add_edge(int u, int v, long long w){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        weight[u][v] = weight[v][u] = max(weight[u][v], w);
    }

    long long maximum_matching(){
        long long max_weight = 0;
        for (auto& row : weight) max_weight = max(max_weight, *max_element(row.begin(), row.end()));

        mate.assign(n, -1);
        dual.assign(2 * n, 0);
        fill(dual.begin(), dual.begin() + n, max_weight);
        top.resize(n);
        iota(top.begin(), top.end(), 0);
        parent.assign(2 * n, -1);
        base.assign(2 * n, -1);
        iota(base.begin(), base.begin() + n, 0);
        label.assign(2 * n, UNLABELED);
        label_edge.assign(2 * n, NONE);
        best_edge.assign(2 * n, NONE);
        children.assign(2 * n, {});
        child_edges.assign(2 * n, {});
        best_edges_to.assign(2 * n, {});
        built_this_stage.assign(2 * n, 0);
        on_path.assign(2 * n, 0);
        unused_ids.clear();
        for (int b = 2 * n - 1; b >= n; b--) unused_ids.push_back(b);

        while (stage()) {}

        long long total = 0;
        for (int v = 0; v < n; v++){
            if (mate[v] > v) total += weight[v][mate[v]];
        }

        return total;
    }

    /// Grows alternating trees from every free vertex, false once the duals prove the matching optimal
    bool stage(){
        fill(label.begin(), label.end(), UNLABELED);
        fill(best_edge.begin(), best_edge.end(), NONE);
        fill(built_this_stage.begin(), built_this_stage.end(), 0);
        queue.clear();
        for (int v = 0; v < n; v++){
            if (mate[v] == -1 && label[top[v]] == UNLABELED) assign_label(v, OUTER, -1);
        }

        while (true){
            for (size_t i = 0; i < queue.size(); i++){
                int v = queue[i];
                for (int u = 0; u < n; u++){
                    if (weight[v][u] <= 0 || top[u] == top[v]) continue;

                    long long s = slack(v, u);
                    int bu = top[u];
                    if (s == 0 && label[bu] == UNLABELED) assign_label(u, INNER, v);
                    else if (s == 0 && label[bu] == OUTER){
                        int b = find_base(v, u);
                        if (b == -1){
                            augment(v, u);
                            end_stage();
                            return true;
                        }
                        add_blossom(b, v, u);
                    }
                    /// Records the tight edge into an inner blossom's vertex, expanding the blossom relabels from it
                    else if (s == 0 && label[u] == UNLABELED) label[u] = INNER, label_edge[u] = {v, u};
                    else if (s > 0 && label[bu] == OUTER){
                        int bv = top[v];
                        if (best_edge[bv] == NONE || s < slack(best_edge[bv])) best_edge[bv] = {v, u};
                    }
                    else if (s > 0 && label[u] == UNLABELED){
                        if (best_edge[u] == NONE || s < slack(best_edge[u])) best_edge[u] = {v, u};
                    }
                }
            }
            queue.clear();

            long long delta = LLONG_MAX;
            int kind = 1, blossom = -1;
            pair<int, int> edge = NONE;
            for (int v = 0; v < n; v++) delta = min(delta, dual[v]);
            for (int v = 0; v < n; v++){
                if (label[top[v]] == UNLABELED && best_edge[v] != NONE && slack(best_edge[v]) < delta){
                    delta = slack(best_edge[v]), kind = 2, edge = best_edge[v];
                }
            }
            for (int b = 0; b < 2 * n; b++){
                if (base[b] == -1 || parent[b] != -1) continue;
                /// Every outer vertex has the parity of the free vertices' duals, so an outer-outer slack is even
                if (label[b] == OUTER && best_edge[b] != NONE && slack(best_edge[b]) / 2 < delta){
                    delta = slack(best_edge[b]) / 2, kind = 3, edge = best_edge[b];
                }
                if (b >= n && label[b] == INNER && dual[b] < delta) delta = dual[b], kind = 4, blossom = b;
            }

            for (int v = 0; v < n; v++){
                if (label[top[v]] == OUTER) dual[v] -= delta;
                else if (label[top[v]] == INNER) dual[v] += delta;
            }
            for (int b = n; b < 2 * n; b++){
                if (base[b] == -1 || parent[b] != -1) continue;
                if (label[b] == OUTER) dual[b] += delta;
                else if (label[b] == INNER) dual[b] -= delta;
            }

            if (kind == 1) return false;
            if (kind == 4) expand(blossom, false);
            else queue.push_back(edge.first);
        }
    }

    void assign_label(int w, int lab, int from){
        int b = top[w];
        label[w] = label[b] = lab;
        label_edge[w] = label_edge[b] = {from, w};
        best_edge[w] = best_edge[b] = NONE;

        if (lab == OUTER){
            for (int x : leaves(b)) queue.push_back(x);
        }
        else assign_label(mate[base[b]], OUTER, base[b]);
    }

    vector<int> leaves(int b){
        vector<int> res, stack = {b};
        while (!stack.empty()){
            int c = stack.back();
            stack.pop_back();
            if (c < n) res.push_back(c);
            else stack.insert(stack.end(), children[c].begin(), children[c].end());
        }

        return res;
    }

    long long slack(int u, int v){
        return dual[u] + dual[v] - 2 * weight[u][v];
    }

    long long slack(pair<int, int> e){
        return slack(e.first, e.second);
    }

    /// Walks up both trees at once: the first blossom seen twice holds the new blossom's base, -1 if the trees differ
    int find_base(int v, int u){
        vector<int> path;
        int res = -1;
        while (v != -1){
            int b = top[v];
            if (on_path[b]){
                res = base[b];
                break;
            }

            on_path[b] = 1;
            path.push_back(b);
            v = label_edge[b].first == -1 ? -1 : label_edge[top[label_edge[b].first]].first;
            if (u != -1) swap(v, u);
        }

        for (int b : path) on_path[b] = 0;
        return res;
    }

    void augment(int v, int u){
        for (auto [s, j] : {pair<int, int>{v, u}, pair<int, int>{u, v}}){
            while (true){
                int bs = top[s];
                if (bs >= n) augment_blossom(bs, s);
                mate[s] = j;
                if (label_edge[bs].first == -1) break;

                int bt = top[label_edge[bs].first];
                tie(s, j) = label_edge[bt];
                if (bt >= n) augment_blossom(bt, j);
                mate[j] = s;
            }
        }
    }

    /// Makes v the base of b by flipping the even path around the cycle from v's child to the old base
    void augment_blossom(int b, int v){
        int t = v;
        while (parent[t] != b) t = parent[t];
        if (t >= n) augment_blossom(t, v);

        auto& ch = children[b];
        auto& ed = child_edges[b];
        int k = ch.size(), i = find(ch.begin(), ch.end(), t) - ch.begin();
        auto match = [&](int j){
            auto [x, y] = ed[j];
            if (ch[j] >= n) augment_blossom(ch[j], x);
            if (ch[(j + 1) % k] >= n) augment_blossom(ch[(j + 1) % k], y);
            mate[x] = y, mate[y] = x;
        };
        if (i % 2 == 1){
            for (int j = i + 1; j < k; j += 2) match(j);
        }
        else{
            for (int j = i - 2; j >= 0; j -= 2) match(j);
        }

        rotate(ch.begin(), ch.begin() + i, ch.end());
        rotate(ed.begin(), ed.begin() + i, ed.end());
        base[b] = v;
    }

    void end_stage(){
        for (int b = n; b < 2 * n; b++){
            if (base[b] != -1 && parent[b] == -1 && label[b] == OUTER && dual[b] == 0) expand(b, true);
        }
    }

    void add_blossom(int base_vertex, int v, int u){
        int bb = top[base_vertex], bv = top[v], bu = top[u];
        int b = unused_ids.back();
        unused_ids.pop_back();
        base[b] = base_vertex, parent[b] = -1, parent[bb] = b;

        vector<int> side_v, side_u;
        vector<pair<int, int>> edges_v, edges_u;
        while (bv != bb){
            parent[bv] = b;
            side_v.push_back(bv);
            edges_v.push_back(label_edge[bv]);
            bv = top[label_edge[bv].first];
        }
        while (bu != bb){
            parent[bu] = b;
            side_u.push_back(bu);
            edges_u.push_back({label_edge[bu].second, label_edge[bu].first});
            bu = top[label_edge[bu].first];
        }

        children[b] = {bb};
        children[b].insert(children[b].end(), side_v.rbegin(), side_v.rend());
        children[b].insert(children[b].end(), side_u.begin(), side_u.end());
        child_edges[b].assign(edges_v.rbegin(), edges_v.rend());
        child_edges[b].push_back({v, u});
        child_edges[b].insert(child_edges[b].end(), edges_u.begin(), edges_u.end());

        label[b] = OUTER, label_edge[b] = label_edge[bb], dual[b] = 0;
        for (int x : leaves(b)){
            if (label[top[x]] == INNER) queue.push_back(x);
            top[x] = b;
        }

        /// One least-slack edge per neighbouring outer blossom, so merging never rescans the leaves of listed children
        vector<pair<int, int>> best_to(2 * n, NONE);
        for (int c : children[b]){
            vector<pair<int, int>> candidates;
            if (built_this_stage[c]) candidates = best_edges_to[c];
            else{
                for (int x : leaves(c)){
                    for (int y = 0; y < n; y++) if (weight[x][y] > 0) candidates.push_back({x, y});
                }
            }

            for (auto [x, y] : candidates){
                int by = top[y];
                if (by != b && label[by] == OUTER && (best_to[by] == NONE || slack(x, y) < slack(best_to[by]))){
                    best_to[by] = {x, y};
                }
            }
        }

        best_edges_to[b].clear();
        best_edge[b] = NONE;
        for (auto e : best_to){
            if (e == NONE) continue;
            best_edges_to[b].push_back(e);
            if (best_edge[b] == NONE || slack(e) < slack(best_edge[b])) best_edge[b] = e;
        }
        built_this_stage[b] = 1;
    }

    /// endstage expands an outer blossom with zero dual between stages, otherwise b is inner and its children relabel
    void expand(int b, bool endstage){
        for (int c : children[b]){
            parent[c] = -1;
            if (c < n) top[c] = c;
            else if (endstage && dual[c] == 0) expand(c, true);
            else{
                for (int x : leaves(c)) top[x] = c;
            }
        }

        if (!endstage && label[b] == INNER){
            auto& ch = children[b];
            auto& ed = child_edges[b];
            int k = ch.size(), entry = top[label_edge[b].second];
            int j = find(ch.begin(), ch.end(), entry) - ch.begin(), dir = j % 2 == 1 ? 1 : k - 1;

            /// The even path from the entry child to the base alternates inner and outer children
            pair<int, int> e = label_edge[b];
            while (j != 0){
                assign_label(e.second, INNER, e.first);
                j = (j + dir) % k;
                int next = (j + dir) % k;
                e = dir == 1 ? ed[j] : pair<int, int>{ed[next].second, ed[next].first};
                j = next;
            }
            label[ch[0]] = label[e.second] = INNER;
            label_edge[ch[0]] = label_edge[e.second] = e;
            best_edge[ch[0]] = NONE;

            /// An odd path child stays unlabeled unless an outer vertex reached one of its vertices by a tight edge
            for (j = (j + dir) % k; ch[j] != entry; j = (j + dir) % k){
                if (label[ch[j]] == OUTER) continue;
                for (int x : leaves(ch[j])){
                    if (label[x] == INNER){
                        assign_label(x, INNER, label_edge[x].first);
                        break;
                    }
                }
            }
        }

        base[b] = -1;
        label[b] = UNLABELED;
        children[b].clear();
        child_edges[b].clear();
        unused_ids.push_back(b);
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

int main(){
    int n, i, j, w;

    read(n);
    WeightedBlossom wb(n);
    for (i = 0; i < n; i++){
        for (j = i + 1; j < n; j++){
            read(w);
            wb.add_edge(i, j, w);
        }
    }

    write(wb.maximum_matching());
    return 0;
}
