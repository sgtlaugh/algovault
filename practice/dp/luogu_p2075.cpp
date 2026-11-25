/***
 *
 * Range LIS
 * Length of the longest increasing subsequence of any subarray, all queries answered in one batch
 *
 * Complexity:
 *   - RangeLIS: O(n log^2 n) to build, O((n + q) log n) for a batch of q queries, O(n) memory
 *   - UnitMongeProduct: O(n log n) per product of two permutations of size n
 *
 * RangeLIS lis(a): a of any type with <, strictly increasing subsequences
 * RangeLIS lis(a, true): non-decreasing subsequences instead
 * lis.query(queries): res[i] = the LIS length of a[l..r] for queries[i] = {l, r}, 0 when l > r
 * UnitMongeProduct product: the permutation product RangeLIS is built on, callable on its own
 * product(a, b): the permutation c with S(c)(i, k) = min over j of S(a)(i, j) + S(b)(j, k)
 *
 * Notes:
 *   - a, b and c are permutations of 0..n-1, and S(p)(i, j) = #{i' >= i : p[i'] < j} for i, j in [0, n]
 *   - n = 10^5 builds in about 0.5 s at -O2
 *   - Source: the seaweed braid of A. Tiskin, Semi-local string comparison (arXiv 0707.3619)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct UnitMongeProduct{
    vector<int> pool;
    size_t top = 0;

    vector<int> operator()(const vector<int>& a, const vector<int>& b){
        int n = a.size();
        assert((int)b.size() == n);
        if (pool.size() < 12 * (size_t)n + 256) pool.resize(12 * (size_t)n + 256);

        vector<int> c(n);
        multiply(a.data(), b.data(), c.data(), n);
        return c;
    }

    /// Splitting the middle index j in two halves, each row keeps its half's point except where the boundary steps
    void multiply(const int* a, const int* b, int* c, int n){
        if (n <= 1){
            if (n) c[0] = 0;
            return;
        }

        int mid = n / 2;
        size_t saved = top;
        int *rows = take(n), *sub_a = take(n), *cols = take(n), *rank = take(n), *sub_b = take(n), *sub_c = take(n);
        int lo = 0, hi = mid;
        for (int i = 0; i < n; i++){
            if (a[i] < mid) rows[lo] = i, sub_a[lo++] = a[i];
            else rows[hi] = i, sub_a[hi++] = a[i] - mid;
        }
        fill(rank, rank + n, 1);
        for (int j = 0; j < mid; j++) rank[b[j]] = 0;
        lo = 0, hi = mid;
        for (int k = 0; k < n; k++){
            if (rank[k] == 0) rank[k] = lo, cols[lo++] = k;
            else rank[k] = hi - mid, cols[hi++] = k;
        }
        for (int j = 0; j < n; j++) sub_b[j] = rank[b[j]];
        multiply(sub_a, sub_b, sub_c, mid);
        multiply(sub_a + mid, sub_b + mid, sub_c + mid, n - mid);

        int *row_of = rank, *in_hi = sub_a, *bnd = take(n + 1);
        for (int t = 0; t < n; t++){
            c[rows[t]] = cols[t < mid ? sub_c[t] : mid + sub_c[t]];
            in_hi[rows[t]] = t >= mid;
        }
        for (int i = 0; i < n; i++) row_of[c[i]] = i;

        /// d(i, k) = #hi points [row < i, col < k] - #lo points [row >= i, col >= k]; bnd[i] = first k with d >= 0
        bnd[n] = 0;
        int d = 0;
        for (int i = n - 1; i >= 0; i--){
            int k = bnd[i + 1];
            d -= in_hi[i] ? c[i] < k : c[i] >= k;
            while (d < 0){
                int r = row_of[k++];
                d += in_hi[r] ? r < i : r >= i;
            }
            bnd[i] = k;
        }
        for (int i = 0; i < n; i++){
            if (bnd[i] > bnd[i + 1]) c[i] = bnd[i] - 1;
        }

        top = saved;
    }

    /// Scratch from one buffer in stack order: a vector per call made the range LIS build 5x slower
    int* take(int len){
        assert(top + len <= pool.size());
        int* res = pool.data() + top;
        top += len;
        return res;
    }
};

struct RangeLIS{
    int n;
    vector<int> start_at;  /// start_at[e]: the top column whose seaweed leaves the bottom at column e, -1 for none
    UnitMongeProduct product;

    template <typename T>
    RangeLIS(const vector<T>& a, bool allow_equal = false) : n(a.size()), start_at(n, -1){
        /// Equal values never rise strictly, and always rise weakly: ranks tie-break by index accordingly
        vector<int> order(n), p(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j){
            if (a[i] < a[j] || a[j] < a[i]) return a[i] < a[j];
            return allow_equal ? i < j : i > j;
        });
        for (int r = 0; r < n; r++) p[order[r]] = r;

        vector<int> out = braid(p);
        for (int s = 0; s < n; s++){
            if (out[n + s] < n) start_at[out[n + s]] = s;
        }
        product.pool = {};
    }

    /// LIS(a[l..r]) = (r - l + 1) - #{seaweeds that start at top column >= l and leave the bottom at column <= r}
    vector<int> query(const vector<pair<int, int>>& queries) const{
        int q = queries.size();
        vector<int> order(q), res(q, 0), tree(n + 1, 0);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j){ return queries[i].second < queries[j].second; });

        int e = 0, added = 0;
        for (int i : order){
            auto [l, r] = queries[i];
            assert(0 <= l && l <= n && -1 <= r && r < n);
            for (; e <= r; e++){
                if (start_at[e] == -1) continue;
                for (int x = start_at[e] + 1; x <= n; x += x & -x) tree[x]++;
                added++;
            }

            int below = 0;
            for (int x = l; x > 0; x -= x & -x) below += tree[x];
            if (l <= r) res[i] = (r - l + 1) - (added - below);
        }

        return res;
    }

    /// p's braid on 2k strands, rows by value: in at left then top, out at bottom then right, each from bottom-left
    vector<int> braid(const vector<int>& p){
        int k = p.size();
        if (k <= 1) return k ? vector<int>{0, 1} : vector<int>{};

        int h = k / 2;
        vector<int> lo_cols, hi_cols, p_lo, p_hi;
        for (int j = 0; j < k; j++){
            if (p[j] < h) lo_cols.push_back(j), p_lo.push_back(p[j]);
            else hi_cols.push_back(j), p_hi.push_back(p[j] - h);
        }
        vector<int> lo = expand(braid(p_lo), lo_cols, k), hi = expand(braid(p_hi), hi_cols, k);

        /// The low rows lie above the high rows; each block passes the other block's left or right strands straight on
        vector<int> above(2 * k), below(2 * k);
        for (int i = 0; i < k - h; i++) above[i] = i;
        for (int i = 0; i < k + h; i++) above[k - h + i] = k - h + lo[i];
        for (int i = 0; i < 2 * k - h; i++) below[i] = hi[i];
        for (int i = 2 * k - h; i < 2 * k; i++) below[i] = i;
        return product(above, below);
    }

    /// A child's braid over the k parent columns: a column the child does not own passes straight down
    static vector<int> expand(const vector<int>& b, const vector<int>& cols, int k){
        int rows = cols.size();
        auto parent_out = [&](int o){ return o < rows ? cols[o] : k + o - rows; };

        vector<int> res(rows + k);
        for (int r = 0; r < rows; r++) res[r] = parent_out(b[r]);
        for (int c = 0; c < k; c++) res[rows + c] = c;
        for (int t = 0; t < rows; t++) res[rows + cols[t]] = parent_out(b[rows + t]);
        return res;
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
    int n, q, i, l, r;

    read(n, q);
    vector<int> ar(n);
    read(ar);

    vector<pair<int, int>> queries(q);
    for (i = 0; i < q; i++){
        read(l, r);
        queries[i] = {l - 1, r - 1};
    }

    RangeLIS lis(ar);
    vector<int> res = lis.query(queries);
    for (i = 0; i < q; i++) write(res[i]);

    return 0;
}
