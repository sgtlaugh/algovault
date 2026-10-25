/***
 *
 * Disjoint Set Union With Rollback
 * Union-find whose unions can be undone, for offline dynamic connectivity and divide and conquer over time
 *
 * Complexity: O(log n) per operation, O(1) per undone union
 *
 * RollbackDSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they are already one set
 * dsu.rollback(t): undo every union made after time() returned t, 0 <= t <= time()
 * dsu.find(a): the root of a's set, the same node for every member
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 * dsu.time(): number of successful connects not yet rolled back, use it as a checkpoint
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED rollback_dsu
struct RollbackDSU{
    vector<int> sz, parent, history;

    RollbackDSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        history.push_back(a);
        return true;
    }

    void rollback(int t){
        assert(0 <= t && t <= time());
        while ((int)history.size() > t){
            int a = history.back();
            history.pop_back();
            sz[parent[a]] -= sz[a], parent[a] = a;
        }
    }

    int find(int i){
        while (i != parent[i]) i = parent[i];
        return i;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }

    int time(){
        return history.size();
    }
};
/// END SHARED rollback_dsu

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

#define MAX 300010
#define OFF 300000

struct state{
    int root, cx, cy;
    long long res;
};

int q, cx[2 * MAX], cy[2 * MAX];
long long res, ans[MAX];
vector<pair<int, int>> tree[4 * MAX];
vector<state> stk;
RollbackDSU dsu(2 * MAX);

void add(int idx, int a, int b, int l, int r, pair<int, int> e){
    if (r < a || b < l) return;
    if (l <= a && b <= r){
        tree[idx].push_back(e);
        return;
    }

    int mid = (a + b) >> 1;
    add(idx * 2, a, mid, l, r, e);
    add(idx * 2 + 1, mid + 1, b, l, r, e);
}

void connect(int u, int v){
    u = dsu.find(u), v = dsu.find(v);
    if (u == v) return;

    dsu.connect(u, v);
    int w = dsu.find(u);
    stk.push_back({w, cx[w], cy[w], res});
    res -= ((long long)cx[u] * cy[u]) + ((long long)cx[v] * cy[v]);
    cx[w] = cx[u] + cx[v], cy[w] = cy[u] + cy[v];
    res += (long long)cx[w] * cy[w];
}

void rollback(int t){
    while ((int)stk.size() > t){
        state s = stk.back();
        stk.pop_back();
        cx[s.root] = s.cx, cy[s.root] = s.cy, res = s.res;
    }
    dsu.rollback(t);
}

void dfs(int idx, int a, int b){
    int t = dsu.time();
    for (auto e: tree[idx]) connect(e.first, e.second);

    if (a == b) ans[a] = res;
    else{
        int mid = (a + b) >> 1;
        dfs(idx * 2, a, mid);
        dfs(idx * 2 + 1, mid + 1, b);
    }
    rollback(t);
}

int main(){
    int i, x, y;

    read(q);
    map<pair<int, int>, int> start;
    for (i = 0; i < q; i++){
        read(x, y);
        auto it = start.find({x, y + OFF});
        if (it == start.end()) start[{x, y + OFF}] = i;
        else{
            add(1, 0, q - 1, it->second, i - 1, it->first);
            start.erase(it);
        }
    }
    for (auto it: start) add(1, 0, q - 1, it.second, q - 1, it.first);

    for (i = 1; i <= OFF; i++) cx[i] = 1, cy[i + OFF] = 1;
    dfs(1, 0, q - 1);

    vector<long long> v(ans, ans + q);
    write(v);
    return 0;
}
