/***
 *
 * Plug DP
 * Counts or optimizes cycles and paths on a grid: the caller writes one cell's transition, the engine runs the rest
 *
 * Complexity: O(S) expected per cell for S states alive at that cell, plus the transition's own cost
 *
 * PlugDP<T> dp(n, m, bits, start): n x m grid, bits per plug, the all-zero state holds start, equal states add up
 * PlugDP<T, Merge> dp(n, m, bits, start, merge): equal states combine their values with merge(a, b), such as min
 * dp.add(s, v): only inside the transition, offers state s with value v after the current cell
 * dp.run(transition): calls transition(row, col, s, v) for each state s of value v before each cell, row-major, once
 * dp.get(s, pos): plug pos of state s
 * dp.match(s, pos): position of the bracket paired with plug pos, plugs 1 opening and 2 closing, balanced
 * dp.normalize(s): s with its nonzero labels renumbered 1, 2, 3... in order of first appearance from position 0
 * dp.result(missing): value of the all-zero state after the last cell, missing (T() by default) if none
 * dp.set(s, pos, x): state s with plug pos replaced by x, 0 <= x < 2^bits
 * dp.states, dp.values: the states alive before the current cell, or after the last one, and their values
 *
 * Notes:
 *   - Before cell (row, col) plug col is its edge from the left and col + 1 from above, after it down and right
 *   - Plugs are edges: add drops states with an edge leaving the grid, so transitions never test the border
 *   - (m + 1) bits <= 64, bits <= 8: 2 bits fit brackets, labels need 2^bits >= components + 2 to leave a fresh one
 *   - A state holds plugs 0..m only: the row shift moves every bit up by bits, so nothing can be kept above plug m
 *   - main() writes Formula 1 (one cycle, bracket encoding) and HDU 1693 (any cycles, 1 bit) as transitions
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Merge = plus<T>>
struct PlugDP{
    int n, m, bits, row = 0, col = 0, lg = 10;
    uint64_t off_grid = 0;
    Merge merge;
    vector<uint64_t> states, next_states;
    vector<T> values, next_values;
    vector<int> head, link;

    PlugDP(int n, int m, int bits, T start, Merge merge = Merge())
        : n(n), m(m), bits(bits), merge(merge), states{0}, values{start}, head(1 << lg, -1){
        assert(n >= 1 && m >= 1 && bits >= 1 && bits <= 8 && (m + 1) * bits <= 64);
    }

    void add(uint64_t s, T v){
        if (s & off_grid) return;

        size_t h = slot(s);
        for (int i = head[h]; i >= 0; i = link[i]){
            if (next_states[i] == s){
                next_values[i] = merge(next_values[i], v);
                return;
            }
        }

        link.push_back(head[h]);
        head[h] = next_states.size();
        next_states.push_back(s);
        next_values.push_back(v);
    }

    template <typename Transition>
    void run(Transition transition){
        uint64_t ones = (uint64_t(1) << bits) - 1;
        for (row = 0; row < n; row++){
            for (col = 0; col < m; col++){
                /// sized per cell, not per insert: a resize check inside add stops GCC inlining it, 1.6x slower
                while (head.size() < 8 * states.size()) head.assign(size_t(1) << ++lg, -1);
                off_grid = (col == m - 1 ? ones << (m * bits) : 0) | (row == n - 1 ? ones << (col * bits) : 0);
                for (size_t i = 0; i < states.size(); i++) transition(row, col, states[i], values[i]);

                for (uint64_t s : next_states) head[slot(s)] = -1;
                link.clear();
                swap(states, next_states), swap(values, next_values);
                next_states.clear(), next_values.clear();
            }
            for (uint64_t& s : states) s <<= bits;
        }
    }

    int get(uint64_t s, int pos) const{
        return s >> (pos * bits) & ((1 << bits) - 1);
    }

    int match(uint64_t s, int pos) const{
        int dir = get(s, pos) == 1 ? 1 : -1, depth = 0;
        for (int p = pos; p >= 0 && p <= m; p += dir){
            int x = get(s, p);
            depth += (x == 1) - (x == 2);
            if (depth == 0) return p;
        }

        assert(false);
        return -1;
    }

    uint64_t normalize(uint64_t s) const{
        unsigned char id[256] = {};
        int labels = 0;
        for (int p = 0; p <= m; p++){
            int x = get(s, p);
            if (!x) continue;
            if (!id[x]) id[x] = ++labels;
            s = set(s, p, id[x]);
        }

        return s;
    }

    T result(T missing = T()) const{
        return states.empty() ? missing : values[0];
    }

    uint64_t set(uint64_t s, int pos, int x) const{
        int shift = pos * bits;
        return (s & ~(uint64_t((1 << bits) - 1) << shift)) | uint64_t(x) << shift;
    }

    size_t slot(uint64_t s) const{
        return s * 0x9E3779B97F4A7C15ULL >> (64 - lg);
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

#define MAX 210

int n, m, md, parent[MAX * 8];
bool right_edge[MAX][8], down_edge[MAX][8];

struct Add{
    long long operator()(long long a, long long b) const{
        return (a + b) % md;
    }
};

int find(int u){
    if (parent[u] == u) return u;

    return (parent[u] = find(parent[u]));
}

bool has_cycle(){
    int i, j, u, v;

    for (i = 0; i < (n * m); i++) parent[i] = i;
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            u = find(i * m + j);
            if (right_edge[i][j]){
                v = find(i * m + j + 1);
                if (u == v) return true;
                parent[v] = u;
            }
            if (down_edge[i][j]){
                v = find((i + 1) * m + j);
                if (u == v) return true;
                parent[v] = u;
            }
        }
    }

    return false;
}

bool has_label(PlugDP<long long, Add>& dp, uint64_t s, int x){
    for (int p = 0; p <= m; p++){
        if (dp.get(s, p) == x) return true;
    }

    return false;
}

uint64_t place(PlugDP<long long, Add>& dp, int j, uint64_t s, int up, int left){
    int p, a = dp.get(s, j + 1), b = j ? dp.get(s, j - 1) : 0, c = 15;
    uint64_t x = dp.set(s, j + 1, 0);

    if ((up && !a) || (left && !b)) return 0;
    if (up && left){
        if (a == b) return 0;
        for (p = 0; p <= m; p++){
            if (dp.get(x, p) == a) x = dp.set(x, p, b);
        }
        c = b;
    }
    else if (up) c = a;
    else if (left) c = b;

    if (!up && a && !has_label(dp, x, a)) return 0;
    return dp.normalize(dp.set(x, j, c));
}

long long solve(){
    PlugDP<long long, Add> dp(n + 1, m, 4, 1 % md);
    vector<unordered_map<uint64_t, array<uint64_t, 4>>> memo(m);

    dp.run([&](int i, int j, uint64_t s, long long v){
        int k, up, left;

        if (i == n){
            if (!has_label(dp, s, 2)) dp.add(0, v);
            return;
        }

        auto it = memo[j].find(s);
        if (it == memo[j].end()){
            array<uint64_t, 4> nxt;
            for (k = 0; k < 4; k++) nxt[k] = place(dp, j, s, k >> 1, k & 1);
            it = memo[j].insert({s, nxt}).first;
        }

        for (up = (i && down_edge[i - 1][j]); up <= (i > 0); up++){
            for (left = (j && right_edge[i][j - 1]); left <= (j > 0); left++){
                k = up * 2 + left;
                if (it->second[k]) dp.add(it->second[k], v);
            }
        }
    });

    return dp.result();
}

int main(){
    int t, i, j;
    string line;

    while (read(t)){
        while (t--){
            read(n, m, md);
            for (i = 0; i < n; i++){
                for (j = 0; j < m; j++) right_edge[i][j] = down_edge[i][j] = false;
            }

            read_line(line);
            if (line.find('.') == string::npos) read_line(line);
            for (i = 0; i < (2 * n - 1); i++){
                if (i) read_line(line);
                for (j = 0; j < (int)line.size(); j++){
                    if ((i & 1) && line[j] == '|') down_edge[i >> 1][j >> 1] = true;
                    if (!(i & 1) && line[j] == '-') right_edge[i >> 1][j >> 1] = true;
                }
            }

            if (has_cycle()) write("Impossible");
            else write(solve());
        }
    }

    return 0;
}
