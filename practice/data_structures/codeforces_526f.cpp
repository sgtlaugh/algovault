/***
 *
 * Permutation Tree
 * Splits a permutation into its common intervals, the ranges whose values form a contiguous set
 *
 * Complexity: O(n log n) to build, O(n) for count_intervals
 *
 * PermutationTree pt(p): the tree of a permutation p of 0..n - 1, n may be 0
 * pt.count_intervals(): number of common intervals, the n singletons included
 * pt.root: the root node, -1 if n = 0
 * pt.parent[v]: parent of v, -1 for the root
 * pt.children[v]: children of v in position order
 * pt.span[v]: the positions [l, r] that v covers
 * pt.range[v]: the values [lo, hi] of p[l..r], always hi - lo == r - l
 * pt.type[v]: INCREASING, DECREASING or CUT
 *
 * Notes:
 *   - Nodes 0..n - 1 are the leaves, node i is position i, internal nodes get ids n, n + 1, ..., at most 2n - 1 in all
 *   - INCREASING and DECREASING are join nodes: children ranges ascend or descend one step each
 *   - CUT is a leaf, or a node with 4+ children where only the whole node is a common interval
 *   - No join node has a child that is a join node of the same direction, which makes the tree unique
 *   - Every common interval is a leaf, a whole cut node, or a run of 2+ consecutive children of a join node
 *   - For 1-indexed or arbitrary distinct values, compress them to ranks 0..n - 1 first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PermutationTree{
    enum Type{CUT, INCREASING, DECREASING};

    int n, root;
    vector<int> parent;
    vector<Type> type;
    vector<pair<int, int>> span, range;
    vector<vector<int>> children;

    PermutationTree(const vector<int>& p) : n(p.size()), root(-1){
        parent.reserve(2 * n), type.reserve(2 * n), span.reserve(2 * n), range.reserve(2 * n), children.reserve(2 * n);
        for (int i = 0; i < n; i++) new_node(CUT, {i, i}, {p[i], p[i]});

        /// tree_min at j holds (max - min of p[j..i]) - (i - j), which is 0 exactly when [j, i] is a common interval
        vector<int> tree_min(4 * n), tree_add(4 * n);
        vector<int> max_stack, min_stack, roots_left_to_right;
        for (int i = 0; i < n; i++){
            while (!max_stack.empty() && p[max_stack.back()] < p[i]){
                int r = max_stack.back();
                max_stack.pop_back();
                add(tree_min, tree_add, 1, 0, n - 1, max_stack.empty() ? 0 : max_stack.back() + 1, r, p[i] - p[r]);
            }
            max_stack.push_back(i);

            while (!min_stack.empty() && p[min_stack.back()] > p[i]){
                int r = min_stack.back();
                min_stack.pop_back();
                add(tree_min, tree_add, 1, 0, n - 1, min_stack.empty() ? 0 : min_stack.back() + 1, r, p[r] - p[i]);
            }
            min_stack.push_back(i);

            int cur = i;
            while (true){
                if (!roots_left_to_right.empty() && (adjacent(roots_left_to_right.back(), cur) || adjacent(cur, roots_left_to_right.back()))){
                    int top = roots_left_to_right.back();
                    Type dir = adjacent(top, cur) ? INCREASING : DECREASING;
                    roots_left_to_right.pop_back();
                    if (type[top] == dir){
                        attach(top, cur);
                        cur = top;
                        continue;
                    }

                    int v = new_node(dir, span[top], range[top]);
                    attach(v, top);
                    attach(v, cur);
                    cur = v;
                    continue;
                }

                int start = span[cur].first;
                if (start == 0 || query_min(tree_min, tree_add, 1, 0, n - 1, 0, start - 1) != 0) break;

                int v = new_node(CUT, span[cur], range[cur]);
                attach(v, cur);
                do {
                    attach(v, roots_left_to_right.back());
                    roots_left_to_right.pop_back();
                } while (range[v].second - range[v].first != span[v].second - span[v].first);
                reverse(children[v].begin(), children[v].end());
                cur = v;
            }
            roots_left_to_right.push_back(cur);

            /// Moving to i + 1 lengthens every [j, i] by one, so the -(i - j) term drops by 1 for every j <= i
            add(tree_min, tree_add, 1, 0, n - 1, 0, i, -1);
        }

        if (n) root = roots_left_to_right.back();
    }

    long long count_intervals() const{
        long long total = 0;
        for (int v = 0; v < (int)children.size(); v++){
            long long k = children[v].size();
            total += type[v] == CUT ? 1 : k * (k - 1) / 2;
        }

        return total;
    }

    int new_node(Type t, pair<int, int> s, pair<int, int> r){
        parent.push_back(-1);
        type.push_back(t);
        span.push_back(s);
        range.push_back(r);
        children.emplace_back();
        return parent.size() - 1;
    }

    static void add(vector<int>& tree_min, vector<int>& tree_add, int idx, int a, int b, int l, int r, int val){
        if (r < a || b < l) return;
        if (l <= a && b <= r){
            tree_min[idx] += val, tree_add[idx] += val;
            return;
        }

        int c = (a + b) >> 1;
        add(tree_min, tree_add, idx << 1, a, c, l, r, val);
        add(tree_min, tree_add, idx << 1 | 1, c + 1, b, l, r, val);
        tree_min[idx] = min(tree_min[idx << 1], tree_min[idx << 1 | 1]) + tree_add[idx];
    }

    /// Order matters: adjacent(left, right) means the merge of left then right is INCREASING
    bool adjacent(int u, int v) const{
        return range[u].second + 1 == range[v].first;
    }

    void attach(int v, int child){
        parent[child] = v;
        children[v].push_back(child);
        span[v] = {min(span[v].first, span[child].first), max(span[v].second, span[child].second)};
        range[v] = {min(range[v].first, range[child].first), max(range[v].second, range[child].second)};
    }

    static int query_min(const vector<int>& tree_min, const vector<int>& tree_add, int idx, int a, int b, int l, int r){
        if (r < a || b < l) return INT_MAX;
        if (l <= a && b <= r) return tree_min[idx];

        int c = (a + b) >> 1;
        return min(query_min(tree_min, tree_add, idx << 1, a, c, l, r),
                   query_min(tree_min, tree_add, idx << 1 | 1, c + 1, b, l, r)) + tree_add[idx];
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
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - Reading into a char array does not compile: read a string, or a grid as vector<string>
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A - not followed by a digit reads as 0 and eats the byte after it, so "3 - 4" reads 3, then 0, then 4
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
 *   - A float, double or long double read gives strtof's, strtod's or strtold's value, 2-3x faster on short decimals
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
        while (isdigit(c)){
            v = v * 10 - (c - '0');

            /// The digits already buffered are read through a local pointer, 15% faster than a read_char call per digit
            const char *p = inbuf + inptr, *end = inbuf + buf_len;
            while (p < end && (unsigned)(*p - '0') < 10) v = v * 10 - (*p++ - '0');
            inptr = p - inbuf;

            /// A buffered non-digit ends the number and is eaten, else the digits may go on after a refill
            if (p < end){
                inptr++;
                break;
            }
            c = read_char();
        }
        x = neg ? v : -v;

        return true;
    }

    /// The C locale's isspace set, false for EOF, and 1.4-1.7x faster for strings and chars than isspace's table lookup
    inline bool is_space(int c){
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (is_space(c)) c = read_char();
        if (c == EOF) return 0;

        for (; !is_space(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return len;
    }

    /// Digits m and 10^|e| both exact in F round once, as parse does: 2-3x faster than parse on short tokens
    template <typename F, F (*parse)(const char *, char **)>
    bool read_floating(F &x, int max_len, unsigned long long max_mantissa, int max_exp){
        int len = read_token();
        if (!len) return false;

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds float and double twice, long double once
        bool own_width = !FLT_EVAL_METHOD || is_same<F, long double>::value;

        /// A longer token skips the digit scan, which would make its parse about 10% slower
        if (len <= max_len && own_width){
            static const F p10[] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L,
                                    1e10L, 1e11L, 1e12L, 1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L,
                                    1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};
            const unsigned char *s = (const unsigned char *)tmpbuf + (*tmpbuf == '-' || *tmpbuf == '+');
            unsigned long long m = 0;
            int digits = 0, e = 0;
            for (; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++;
            if (*s == '.'){
                for (s++; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++, e--;
            }

            /// No exponent digits after the e equals strtod, which then ignores the e
            if (*s == 'e' || *s == 'E'){
                int sign = s[1] == '-' ? -1 : 1, ex = 0;
                s += 1 + (s[1] == '-' || s[1] == '+');
                for (; isdigit(*s) && ex < 1000; s++) ex = ex * 10 + (*s - '0');
                e += sign * ex;
            }

            if (!*s && digits && digits <= 19 && m <= max_mantissa && e >= -max_exp && e <= max_exp){
                F v = e < 0 ? m / p10[-e] : m * p10[e];
                x = *tmpbuf == '-' ? -v : v;
                return true;
            }
        }

        x = parse(tmpbuf, nullptr);
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        return read_floating<float, strtof>(x, 10, 1 << 24, 10);
    }

    /// 22 characters keep 16-digit scientific tokens like -1.234567890123456e-05 on the fast path
    bool read_one(double &x){
        return read_floating<double, strtod>(x, 22, 1ULL << 53, 22);
    }

    /// Only an x87 64-bit mantissa takes the fast path: a double-double long double does not round a division once
    bool read_one(long double &x){
        return read_floating<long double, strtold>(x, 21, ~0ULL, LDBL_MANT_DIG == 64 ? 27 : -1);
    }

    bool read_one(string &s){
        int c = ' ';

        while (is_space(c)) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !is_space(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (is_space(c));
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

    /// Without it a char * picks the integer template, an exact match that beats the const char * overload
    void write_one(char *s){
        write_one((const char *)s);
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
        static const long long p10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
            10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 1000000000000000,
            10000000000000000, 100000000000000000, 1000000000000000000};

        if (f.digits >= 0 && f.digits <= 18){
            long double scaled = fabsl(f.value) * p10[f.digits];

            /// Rounding never moves the product across k + 0.5, exact below 2^52, so only an exact k + 0.5 needs printf
            if (scaled < 4503599627370496.0L && scaled - (long long)scaled != 0.5L){
                long long q = scaled + 0.5L;
                if (signbit(f.value)) write_char('-');
                write_one(q / p10[f.digits]);
                if (!f.digits) return;

                long long r = q % p10[f.digits];
                for (int j = f.digits - 1; j >= 0; j--) tmpbuf[j] = '0' + r % 10, r /= 10;
                write_char('.');
                for (int j = 0; j < f.digits; j++) write_char(tmpbuf[j]);
                return;
            }
        }

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
    int n, i, r, c;

    read(n);
    vector<int> p(n);
    for (i = 0; i < n; i++){
        read(r, c);
        p[r - 1] = c - 1;
    }

    PermutationTree pt(p);
    write(pt.count_intervals());
    return 0;
}
