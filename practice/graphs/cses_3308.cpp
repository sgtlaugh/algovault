/***
 *
 * Chromatic Number
 * Minimum number of colors in a proper vertex coloring of an undirected graph
 *
 * Complexity: O(2^n n) time, 12 * 2^n bytes of memory (about 200 MB at n = 24)
 *
 * ChromaticNumber g(n): vertices 0..n - 1, n <= 24
 * g.add_edge(u, v): undirected edge, u != v, duplicates are fine
 * g.solve(): the chromatic number, 0 for the empty graph
 *
 * Notes:
 *   - Counts the k-tuples of independent sets covering every vertex mod the prime 2^61 - 1, so it is never too small
 *   - It is too large only if that count at the true chromatic number is a multiple of 2^61 - 1
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ChromaticNumber{
    static constexpr uint64_t MOD = (1ULL << 61) - 1;

    int n;
    vector<int> adj;

    ChromaticNumber(int n) : n(n), adj(n, 0){
        assert(0 <= n && n <= 24);
    }

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        adj[u] |= 1 << v;
        adj[v] |= 1 << u;
    }

    int solve() const{
        if (n == 0) return 0;
        int full = 1 << n;

        vector<uint32_t> ind(full);
        ind[0] = 1;
        for (int mask = 1; mask < full; mask++){
            int u = __builtin_ctz(mask), rest = mask ^ (1 << u);
            ind[mask] = ind[rest] + ind[rest & ~adj[u]];
        }

        vector<uint64_t> power(full, 1);
        for (int k = 1; k < n; k++){
            uint64_t by_parity[2] = {0, 0};
            for (int mask = 0; mask < full; mask++){
                power[mask] = mul(power[mask], ind[mask]);
                uint64_t& sum = by_parity[__builtin_popcount(mask) & 1];
                sum += power[mask];
                if (sum >= MOD) sum -= MOD;
            }
            if (by_parity[0] != by_parity[1]) return k;
        }
        return n;
    }

    static uint64_t mul(uint64_t a, uint64_t b){
        unsigned __int128 product = (unsigned __int128)a * b;
        uint64_t res = (uint64_t)(product & MOD) + (uint64_t)(product >> 61);
        return res >= MOD ? res - MOD : res;
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

int n, adj[16], color[16];
vector<int> cls;

int chromatic(int from){
    vector<int> nodes = cls;

    for (int i = from; i < n; i++) nodes.push_back(1 << i);

    int i, j, k = nodes.size();
    ChromaticNumber g(k);
    for (i = 0; i < k; i++){
        int nb = 0;
        for (j = 0; j < n; j++){
            if ((nodes[i] >> j) & 1) nb |= adj[j];
        }
        for (j = i + 1; j < k; j++){
            if ((nb & nodes[j]) || j < (int)cls.size()) g.add_edge(i, j);
        }
    }

    return g.solve();
}

int main(){
    int m, i, j, a, b, k;

    read(n, m);
    for (i = 0; i < m; i++){
        read(a, b);
        a--, b--;
        adj[a] |= 1 << b, adj[b] |= 1 << a;
    }

    k = chromatic(0);
    for (i = 0; i < n; i++){
        for (j = 0; j < (int)cls.size(); j++){
            if (cls[j] & adj[i]) continue;
            cls[j] |= 1 << i;
            if (chromatic(i + 1) == k) break;
            cls[j] ^= 1 << i;
        }
        if (j == (int)cls.size()) cls.push_back(1 << i);
        color[i] = j + 1;
    }

    write(k);
    write(vector<int>(color, color + n));
    return 0;
}
