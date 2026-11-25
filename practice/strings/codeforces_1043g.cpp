/***
 *
 * Main-Lorentz Tandem Repeats
 * Finds every tandem repeat (square) s[i..i + 2l) with s[i..i + l) == s[i + l..i + 2l), grouped into compact triples
 *
 * Complexity: O(n log n) time and at most n * ceil(log2 n) triples, O(n) extra memory
 *
 * tandem_repeats(s, report): calls report(first, last, l) per triple, a square s[i..i + 2l) at every first <= i <= last
 *
 * Notes:
 *   - Every tandem repeat belongs to exactly one triple, triples come in no particular order
 *   - The repeats can number Theta(n^2) ("aaaa..."): count them as the sum of last - first + 1 in a long long
 *   - Works on strings, vectors and any indexable container whose elements compare with ==, of any values
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace main_lorentz{
    /// z[i] = longest common prefix of s(i..) and s(0..), s maps an index to an element
    template <typename S>
    void z_function(S s, int n, int* z){
        if (n) z[0] = n;

        for (int i = 1, l = 0, r = 0; i < n; i++){
            int val = i < r ? min(r - i, z[i - l]) : 0;
            while (i + val < n && s(val) == s(i + val)) val++;
            if (i + val > r) l = i, r = i + val;
            z[i] = val;
        }
    }

    /// res[i] = longest common prefix of text(i..) and pattern(0..), zp = z of pattern
    template <typename P, typename S>
    void prefix_matches(P pattern, int m, const int* zp, S text, int n, int* res){
        for (int i = 0, l = 0, r = 0; i < n; i++){
            int val = i < r ? min(r - i, zp[i - l]) : 0;
            while (i + val < n && val < m && text(i + val) == pattern(val)) val++;
            if (i + val > r) l = i, r = i + val;
            res[i] = val;
        }
    }

    /// Reports the repeats inside s[lo..hi), buf holds 2 * (hi - lo) ints the children reuse as they finish first
    template <typename Container, typename F>
    void solve(const Container& s, int lo, int hi, int* buf, F& report){
        if (hi - lo < 2) return;

        int mid = lo + (hi - lo) / 2, n_left = mid - lo, n_right = hi - mid;
        solve(s, lo, mid, buf, report);
        solve(s, mid, hi, buf, report);

        /// decltype(auto) keeps const T& for ordinary containers and returns vector<bool>'s proxy by value
        auto left = [&](int i) -> decltype(auto){ return s[lo + i]; };
        auto right = [&](int i) -> decltype(auto){ return s[mid + i]; };
        auto left_reversed = [&](int i) -> decltype(auto){ return s[mid - 1 - i]; };
        auto right_reversed = [&](int i) -> decltype(auto){ return s[hi - 1 - i]; };
        int *z_left_reversed = buf, *left_vs_right = buf + n_left;
        int *z_right = buf + 2 * n_left, *right_reversed_vs_left_reversed = buf + 2 * n_left + n_right;
        z_function(left_reversed, n_left, z_left_reversed);
        z_function(right, n_right, z_right);
        prefix_matches(right, n_right, z_right, left, n_left, left_vs_right);
        prefix_matches(left_reversed, n_left, z_left_reversed, right_reversed, n_right, right_reversed_vs_left_reversed);

        /// Second half starts in left: the repeat holds the pair (mid - l, mid) and l1 >= 1 matches before mid - l
        for (int l = 1; l <= n_left; l++){
            int back = l < n_left ? z_left_reversed[l] : 0, front = left_vs_right[n_left - l];
            int lo1 = max(1, l - front), hi1 = min(l - 1, back);
            if (lo1 <= hi1) report(mid - l - hi1, mid - l - lo1, l);
        }

        /// Second half starts in right: the repeat holds pair (mid - 1, mid - 1 + l), l1 >= 1 matches up to mid - 1
        for (int l = 1; l <= n_right; l++){
            int back = right_reversed_vs_left_reversed[n_right - l], front = l < n_right ? z_right[l] : 0;
            int lo1 = max(1, l - front), hi1 = min(l, back);
            if (lo1 <= hi1) report(mid - hi1, mid - lo1, l);
        }
    }
}

template <typename Container, typename F>
void tandem_repeats(const Container& s, F report){
    int n = s.size();
    vector<int> buf(2 * n);
    main_lorentz::solve(s, 0, n, buf.data(), report);
}

template <typename F>
void tandem_repeats(const char* s, F report){
    tandem_repeats(string(s), report);
}

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

#define MAX 200010
#define LOG 18

const int B = 640;
const int INF = 1e9;
const unsigned long long HMOD = (1ULL << 61) - 1;

string str;
int n, sq_start[MAX], sq_end[MAX], nxt[MAX], prv[MAX], spf[MAX], parent[MAX];
int rmq[LOG][MAX];
unsigned long long base, pw[MAX], h[MAX];
vector<pair<int, int>> by_len[MAX];
vector<pair<unsigned long long, int>> windows;

unsigned long long mul(unsigned long long a, unsigned long long b){
    unsigned __int128 c = (unsigned __int128)a * b;
    unsigned long long x = (unsigned long long)(c & HMOD) + (unsigned long long)(c >> 61);
    return x >= HMOD ? x - HMOD : x;
}

unsigned long long get_hash(int l, int len){
    unsigned long long x = h[l + len], y = mul(h[l], pw[len]);
    return x >= y ? x - y : x + HMOD - y;
}

bool is_equal(int a, int b, int len){
    return get_hash(a, len) == get_hash(b, len);
}

int find(int x){
    while (parent[x] != x) parent[x] = parent[parent[x]], x = parent[x];
    return x;
}

void paint(int* val, bool at_end){
    int i, len, x;

    for (i = 0; i <= n; i++) parent[i] = i, val[i] = INF;
    for (len = 1; (len * 2) <= n; len++){
        int d = at_end ? (2 * len - 1) : 0;
        for (auto p : by_len[len]){
            for (x = find(p.first + d); x <= (p.second + d); x = find(x)){
                val[x] = len;
                parent[x] = x + 1;
            }
        }
    }
}

void gen(){
    int i, j;

    for (i = 2; i <= n; i++){
        if (spf[i]) continue;
        for (j = i; j <= n; j += i){
            if (!spf[j]) spf[j] = i;
        }
    }

    base = mt19937_64(chrono::steady_clock::now().time_since_epoch().count())() % (HMOD - 1000) + 500;
    pw[0] = 1;
    for (i = 0; i < n; i++){
        pw[i + 1] = mul(pw[i], base);
        h[i + 1] = (mul(h[i], base) + (str[i] - 'a' + 1)) % HMOD;
    }
    for (i = 0; (i + B) <= n; i++) windows.push_back({get_hash(i, B), i});
    sort(windows.begin(), windows.end());

    vector<int> last(26, n);
    for (i = n - 1; i >= 0; i--){
        nxt[i] = last[str[i] - 'a'];
        last[str[i] - 'a'] = i;
    }
    last.assign(26, -1);
    for (i = 0; i < n; i++){
        prv[i] = last[str[i] - 'a'];
        last[str[i] - 'a'] = i;
    }

    tandem_repeats(str, [&](int first, int last, int len){ by_len[len].push_back({first, last}); });
    paint(sq_start, false);
    paint(sq_end, true);

    for (i = 0; i < n; i++) rmq[0][i] = sq_start[i] == INF ? INF : (i + 2 * sq_start[i] - 1);
    for (j = 1; j < LOG; j++){
        for (i = 0; (i + (1 << j)) <= n; i++) rmq[j][i] = min(rmq[j - 1][i], rmq[j - 1][i + (1 << (j - 1))]);
    }
}

int query_min(int l, int r){
    int k = __lg(r - l + 1);
    return min(rmq[k][l], rmq[k][r - (1 << k) + 1]);
}

bool has_repeat(int l, int r){
    if ((r - l + 1) > 26) return true;

    vector<bool> seen(26, false);
    for (int i = l; i <= r; i++){
        if (seen[str[i] - 'a']) return true;
        seen[str[i] - 'a'] = true;
    }

    return false;
}

bool is_power(int l, int len){
    for (int x = len; x > 1; ){
        int p = spf[x];
        if (is_equal(l, l + len / p, len - len / p)) return true;
        while ((x % p) == 0) x /= p;
    }

    return false;
}

bool has_border(int l, int r){
    int len = r - l + 1;

    for (int i = 1; i < B && i <= (len >> 1); i++){
        if (is_equal(l, r - i + 1, i)) return true;
    }
    if (len < (2 * B)) return false;

    unsigned long long x = get_hash(l, B);
    auto it = lower_bound(windows.begin(), windows.end(), make_pair(x, r - (len >> 1) + 1));
    for (; it != windows.end() && it->first == x && it->second <= (r - B + 1); it++){
        if (is_equal(l, it->second, r - it->second + 1)) return true;
    }

    return false;
}

int solve(int l, int r){
    int len = r - l + 1;

    if (!has_repeat(l, r)) return -1;
    if (is_power(l, len)) return 1;
    if ((2 * sq_start[l]) <= len || (2 * sq_end[r]) <= len || has_border(l, r)) return 2;
    if (nxt[l] <= r || prv[r] >= l || query_min(l, r) <= r) return 3;
    return 4;
}

int main(){
    int q, l, r;

    read(n, str, q);
    gen();

    while (q--){
        read(l, r);
        write(solve(l - 1, r - 1));
    }

    return 0;
}
