/***
 *
 * 2SAT with Kosaraju's algorithm
 * Finds a true or false value for every variable that satisfies a system of constraints on pairs of variables
 *
 * Complexity: O(n + m), n = number of variables, m = number of clauses
 *
 * TwoSat g(n): no constraints on variables 1 to n, literal x means x is true and -x means x is false
 * g.add_and(a, b): requires both literals a and b
 * g.add_implication(a, b): requires b whenever a holds
 * g.add_or(a, b): requires at least one of a and b
 * g.add_var(): adds variable n + 1 and returns its index
 * g.add_xor(a, b): requires exactly one of a and b
 * g.at_most_one(literals): allows at most one of the literals to hold, a literal listed twice counts twice
 * g.force_false(x): requires literal x to be false
 * g.force_true(x): requires literal x to be true
 * g.get_assignment(): the true variables, ascending, valid after is_satisfiable() returns true
 * g.is_satisfiable(): whether the constraints can all hold, finding an assignment when they can
 * g.value(x): whether literal x holds in that assignment, valid after is_satisfiable() returns true
 *
 * Notes:
 *   - at_most_one adds k - 2 variables via add_var() and 3k - 5 clauses for k >= 2, not k^2 / 2 pairwise clauses
 *   - Variables from add_var() are numbered after the existing ones, so earlier indices never change
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TwoSat{
    int n, t;
    vector<pair<int, int>> edges;  /// signed literals, mapped to nodes 1..2n only in build() since add_var() changes n
    vector<int> adj_start, adj_to, rev_start, rev_to;
    vector<int> visited, comp, order, dfs_t;

    TwoSat(int n): n(n) {}

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    inline void add_implication(int a, int b){
        assert(a != 0 && abs(a) <= n && b != 0 && abs(b) <= n);
        edges.push_back({a, b}), edges.push_back({-b, -a});
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    int add_var(){
        return ++n;
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    /// prev is a literal meaning "one of lits[0..i-1] is true", it forbids lits[i] and carries over to the next prefix
    void at_most_one(const vector<int>& lits){
        int k = lits.size();
        if (k <= 1) return;

        int prev = lits[0];
        for (int i = 1; i + 1 < k; i++){
            int cur = add_var();
            add_implication(prev, -lits[i]);
            add_implication(prev, cur);
            add_implication(lits[i], cur);
            prev = cur;
        }

        add_implication(prev, -lits[k - 1]);
    }

    inline void force_false(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({x, -x});  /// x -> !x is its own contrapositive, add it once
    }

    inline void force_true(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({-x, x});  /// !x -> x is its own contrapositive, add it once
    }

    vector<int> get_assignment(){
        vector<int> res;
        for (int x = 1; x <= n; x++){
            if (value(x)) res.push_back(x);
        }

        return res;
    }

    bool is_satisfiable(){
        build();
        for (int i = 1; i <= n; i++){
            if (comp[i] == comp[i + n]) return false;
        }

        return true;
    }

    /// If x implies !x then comp[!x] is discovered first, so x is false
    inline bool value(int x){
        if (x < 0) return !value(-x);
        return comp[x] < comp[x + n];
    }

    /// Components are numbered in discovery order, which is a reverse topological order of the implication graph
    void build(){
        int m = 2 * n + 2, c = 0;
        flatten(adj_start, adj_to, false);
        flatten(rev_start, rev_to, true);
        dfs_t.assign(m, 0), order.assign(m, 0), comp.assign(m, 0);

        vector<pair<int, int>> topsort_stack;
        visited.assign(m, 0);
        t = 0;
        for (int i = 2 * n; i >= 1; i--){
            if (!visited[i]) topsort(i, topsort_stack);
            order[dfs_t[i]] = i;
        }

        vector<int> dfs_stack;
        visited.assign(m, 0);
        for (int i = 2 * n; i >= 1; i--){
            int x = order[i];
            if (!visited[x]) dfs(x, c++, dfs_stack);
        }
    }

    void flatten(vector<int>& start, vector<int>& to, bool reverse){
        start.assign(2 * n + 2, 0), to.resize(edges.size());
        for (auto [a, b]: edges) start[node(reverse ? b : a)]++;
        for (int i = 1; i <= 2 * n + 1; i++) start[i] += start[i - 1];
        for (auto [a, b]: edges){
            int u = node(a), v = node(b);
            if (reverse) swap(u, v);
            to[--start[u]] = v;
        }
    }

    /// Literal x is node x, literal -x is node n + x
    inline int node(int x){
        return x > 0 ? x : n - x;
    }

    void topsort(int s, vector<pair<int, int>>& stack){
        visited[s] = true, stack.push_back({s, rev_start[s]});
        while (!stack.empty()){
            auto& [u, e] = stack.back();
            if (e < rev_start[u + 1]){
                int v = rev_to[e++];
                if (!visited[v]) visited[v] = true, stack.push_back({v, rev_start[v]});
            }
            else{
                dfs_t[u] = ++t;
                stack.pop_back();
            }
        }
    }

    /// Labels everything reachable from s, the visiting order does not matter here
    void dfs(int s, int c, vector<int>& stack){
        comp[s] = c, visited[s] = true, stack.push_back(s);
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            for (int e = adj_start[u]; e < adj_start[u + 1]; e++){
                int v = adj_to[e];
                if (!visited[v]) comp[v] = c, visited[v] = true, stack.push_back(v);
            }
        }
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

#define MAXV 2000010

int spf[MAXV];

void sieve(){
    for (int i = 2; i < MAXV; i++){
        if (spf[i]) continue;
        for (int j = i; j < MAXV; j += i){
            if (!spf[j]) spf[j] = i;
        }
    }
}

void add(vector<pair<int, int>>& v, int x, int lit){
    while (x > 1){
        int p = spf[x];
        v.push_back({p, lit});
        while (x % p == 0) x /= p;
    }
}

int main(){
    int n, i, j, k, a, b;

    sieve();
    read(n);
    vector<pair<int, int>> v;
    for (i = 1; i <= n; i++){
        read(a, b);
        add(v, a, -i);
        add(v, b, i);
    }

    sort(v.begin(), v.end());
    TwoSat g(n);
    for (i = 0; i < (int)v.size(); i = j){
        vector<int> lits;
        for (j = i; j < (int)v.size() && v[j].first == v[i].first; j++) lits.push_back(v[j].second);
        g.at_most_one(lits);
    }

    k = g.is_satisfiable();
    write(k ? "Yes" : "No");
    return 0;
}
