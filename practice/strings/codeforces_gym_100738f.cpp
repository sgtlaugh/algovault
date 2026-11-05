/***
 *
 * Suffix Array (SA-IS)
 * Sorted suffixes and their LCP array in O(n), with O(1) common prefix and pattern occurrence queries
 *
 * Complexity:
 *   - suffix_array: O(n) when the values span at most n + 256, else O(n log n) to compress them
 *   - SuffixArrayQueries: O(n log n) more for a sparse table, about 76 MB at n = 1e6, plus a copy of the text
 *   - common_prefix: O(1)
 *   - occurrences(p, len): O(log n)
 *   - occurrences(pattern): O(|pattern| log n)
 *
 * SuffixArray res = suffix_array(c): c is a string or any container of integers, of any values
 * res.sa[i]: start of the i-th smallest suffix
 * res.rank[p]: position of the suffix starting at p in sa
 * res.lcp[i]: longest common prefix of the suffixes at sa[i] and sa[i + 1], lcp[n - 1] = 0
 * res.distinct_substrings(): number of distinct non-empty substrings, O(n)
 *
 * SuffixArrayQueries q(c): everything in SuffixArray plus the calls below
 * q.common_prefix(i, j): longest common prefix of the suffixes starting at i and j
 * q.occurrences(p, len): [lo, hi) with sa[lo .. hi) the starts of s[p .. p + len), 0 <= p, 0 <= len, p + len <= n
 * q.occurrences(pattern): the same range for any pattern
 *
 * Notes:
 *   - The count is hi - lo, an absent pattern gives lo == hi at its insertion point, an empty one gives [0, n)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace SuffixArraySAIS{
    /// Sorts the suffixes of text, values in [0, alphabet), as if a sentinel below every value ended it
    template<typename Char>
    void sais(const Char* text, int* sa, int n, int alphabet){
        if (n == 1){
            sa[0] = 0;
            return;
        }

        vector<int> start(alphabet + 1, 0), bucket(alphabet);  /// bucket c spans sa[start[c] .. start[c + 1])
        for (int i = 0; i < n; i++) start[text[i] + 1]++;
        for (int c = 0; c < alphabet; c++) start[c + 1] += start[c];

        /// LMS: an S-type suffix (smaller than the next one) right after an L-type one, visited right to left
        auto for_each_lms = [&](auto visit){
            bool next_s = false;
            for (int i = n - 2; i >= 0; i--){
                bool s = text[i] < text[i + 1] || (text[i] == text[i + 1] && next_s);
                if (!s && next_s) visit(i + 1);
                next_s = s;
            }
        };

        /// ~p marks p - 1 for the other pass, no type array: p - 1 has the type of p when text[p - 1] == text[p]
        auto induce = [&](){
            copy(start.begin(), start.end() - 1, bucket.begin());
            auto place_l = [&](int p){ sa[bucket[text[p]]++] = p > 0 && text[p - 1] < text[p] ? ~p : p; };
            place_l(n - 1);
            for (int i = 0; i < n; i++){
                if (sa[i] > 0) place_l(sa[i] - 1);
            }

            copy(start.begin() + 1, start.end(), bucket.begin());
            for (int i = n - 1; i >= 0; i--){
                if (sa[i] >= 0) continue;

                sa[i] = ~sa[i];
                int p = sa[i] - 1;
                sa[--bucket[text[p]]] = p > 0 && text[p - 1] <= text[p] ? ~p : p;
            }
        };

        /// Inducing from unsorted LMS seeds sorts the LMS substrings
        fill(sa, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        int lms_count = 0;
        for_each_lms([&](int p){
            sa[--bucket[text[p]]] = p;
            lms_count++;
        });
        induce();
        if (lms_count == 0) return;  /// non-increasing text, the L pass sorted it

        /// Compacts the LMS positions to the front of sa, in their sorted order
        int sorted = 0;
        for (int i = 0; i < n; i++){
            int p = sa[i];
            if (p == 0 || text[p - 1] <= text[p]) continue;

            int end = p + 1;
            while (end < n && text[end] == text[p]) end++;
            if (end < n && text[end] > text[p]) sa[sorted++] = p;
        }

        /// LMS positions are at least 2 apart, so sa[lms_count + p / 2] is a free slot for its length, then its name
        fill(sa + lms_count, sa + n, 0);
        int next = 0;
        for_each_lms([&](int p){
            sa[lms_count + p / 2] = next ? next - p : 0;  /// the last one ends at the sentinel and equals no other
            next = p;
        });

        /// Names start at 1, as 0 marks a free slot
        int names = 0, prev = -1, prev_len = 0;
        for (int i = 0; i < lms_count; i++){
            int p = sa[i], len = sa[lms_count + p / 2];
            bool same = prev >= 0 && len == prev_len;
            for (int d = 0; same && d <= len; d++) same = text[p + d] == text[prev + d];
            if (!same) names++;

            sa[lms_count + p / 2] = names;
            prev = p, prev_len = len;
        }

        /// Names in text order sort like the LMS suffixes, lms_count <= n / 2 so they and their suffix array fit in sa
        int* reduced = sa + n - lms_count;
        int filled = lms_count;
        for (int i = n - 1; i >= lms_count; i--){
            if (sa[i] > 0) reduced[--filled] = sa[i] - 1;
        }

        if (names == lms_count){
            for (int i = 0; i < lms_count; i++) sa[reduced[i]] = i;
        }
        else{
            sais(reduced, sa, lms_count, names);
        }

        filled = lms_count;
        for_each_lms([&](int p){ reduced[--filled] = p; });
        for (int i = 0; i < lms_count; i++) sa[i] = reduced[sa[i]];

        fill(sa + lms_count, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        for (int i = lms_count - 1; i >= 0; i--){  /// the i-th smallest lands at or after i, past every unread seed
            int p = sa[i];
            sa[i] = 0;
            sa[--bucket[text[p]]] = p;
        }
        induce();
    }

    template<typename Container>
    vector<int> build_lcp(const Container& c, const vector<int>& sa, const vector<int>& rank){
        int n = sa.size(), k = 0;
        vector<int> lcp(n);

        for (int i = 0; i < n; i++){
            if (rank[i] == n - 1){
                k = 0;
                continue;
            }

            int j = sa[rank[i] + 1];
            while (i + k < n && j + k < n && c[i + k] == c[j + k]) k++;
            lcp[rank[i]] = k;
            if (k > 0) k--;  /// Kasai: dropping the first character loses at most 1 of the match, so O(n) total
        }

        return lcp;
    }
}

struct SuffixArray{
    int n;
    vector<int> sa, rank, lcp;

    long long distinct_substrings() const {
        long long total = ((long long)n * (n + 1)) / 2;
        for (int i = 0; i + 1 < n; i++) total -= lcp[i];
        return total;
    }
};

template<typename Container>
SuffixArray suffix_array(const Container& c){
    int n = c.size();
    SuffixArray result;

    result.n = n;

    if (n == 0) return result;

    /// Values are ranked in T, a cast to long long would wrap unsigned values past LLONG_MAX
    using T = typename Container::value_type;
    using U = make_unsigned_t<T>;
    T lo = *min_element(c.begin(), c.end()), hi = *max_element(c.begin(), c.end());
    auto offset = [&](T x){ return (unsigned long long)U(U(x) - U(lo)); };  /// exact as x >= lo

    bool compress = offset(hi) > (unsigned long long)n + 256;
    vector<T> values;
    if (compress){
        values.assign(c.begin(), c.end());
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());
    }
    auto value_rank = [&](T x){ return compress ? int(lower_bound(values.begin(), values.end(), x) - values.begin()) : int(offset(x)); };
    int alphabet = compress ? values.size() : offset(hi) + 1;

    result.sa.resize(n);
    auto build = [&](auto text){
        for (int i = 0; i < n; i++) text[i] = value_rank(c[i]);
        SuffixArraySAIS::sais(text.data(), result.sa.data(), n, alphabet);
    };

    /// Byte text measured 30% faster than int at n = 1e7, the induce passes read it at random
    if (alphabet <= 256) build(vector<unsigned char>(n));
    else build(vector<int>(n));

    result.rank.resize(n);
    for (int i = 0; i < n; i++) result.rank[result.sa[i]] = i;
    result.lcp = SuffixArraySAIS::build_lcp(c, result.sa, result.rank);

    return result;
}

SuffixArray suffix_array(const char* s){
    return suffix_array(string(s));
}

template<typename Container>
struct SuffixArrayQueries : SuffixArray{
    Container text;
    vector<vector<int>> table;  /// table[k][i] = min(lcp[i .. i + 2^k))

    SuffixArrayQueries(const Container& c) : SuffixArray(suffix_array(c)), text(c){
        table.push_back(lcp);
        for (int k = 1; (1 << k) <= n; k++){
            vector<int> level(n - (1 << k) + 1);
            for (int i = 0; i < (int)level.size(); i++) level[i] = min(table[k - 1][i], table[k - 1][i + (1 << (k - 1))]);
            table.push_back(move(level));
        }
    }

    int common_prefix(int i, int j) const {
        if (i == j) return n - i;

        int a = rank[i], b = rank[j];
        if (a > b) swap(a, b);
        return lcp_min(a, b);
    }

    pair<int, int> occurrences(int p, int len) const {
        if (len == 0) return {0, n};

        int r = rank[p], lo = 0, hi = r;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (lcp_min(mid, r) >= len) hi = mid;
            else lo = mid + 1;
        }

        int first = lo;
        lo = r, hi = n - 1;
        while (lo < hi){
            int mid = (lo + hi + 1) / 2;
            if (lcp_min(r, mid) >= len) lo = mid;
            else hi = mid - 1;
        }

        return {first, lo + 1};
    }

    pair<int, int> occurrences(const Container& pattern) const {
        return {bound(pattern, false), bound(pattern, true)};
    }

    int lcp_min(int l, int r) const {
        int k = __lg(r - l);
        return min(table[k][l], table[k][r - (1 << k)]);
    }

    /// First rank whose suffix cut to |pattern| is >= pattern, or > when strict
    int bound(const Container& pattern, bool strict) const {
        int lo = 0, hi = n;
        while (lo < hi){
            int mid = (lo + hi) / 2, cmp = compare(sa[mid], pattern);
            if (cmp < 0 || (strict && cmp == 0)) lo = mid + 1;
            else hi = mid;
        }

        return lo;
    }

    /// Sign of the suffix at p cut to |pattern| minus pattern, a suffix that ends first is smaller
    int compare(int p, const Container& pattern) const {
        int m = pattern.size();
        for (int k = 0; k < m; k++){
            if (p + k == n) return -1;
            if (text[p + k] != pattern[k]) return text[p + k] < pattern[k] ? -1 : 1;
        }

        return 0;
    }
};

SuffixArrayQueries(const char*) -> SuffixArrayQueries<string>;

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
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
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

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
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

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
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

#define MAX 100010

int n, q, len[MAX], k[MAX], order[MAX], res[MAX], bit[MAX];

void add(int i){
    for (i++; i <= n; i += (i & -i)) bit[i]++;
}

int kth(int x){
    int i, pos = 0;

    for (i = 1 << 17; i; i >>= 1){
        if ((pos + i) <= n && bit[pos + i] < x){
            pos += i;
            x -= bit[pos];
        }
    }

    return pos;
}

int main(){
    int i, j;
    string s;

    read_line(s);
    n = s.size();
    read(q);
    for (i = 0; i < q; i++){
        read(len[i], k[i]);
        order[i] = i;
    }

    SuffixArray sa = suffix_array(s);
    sort(order, order + q, [](int a, int b){ return len[a] > len[b]; });

    for (i = 0, j = 0; i < q; i++){
        int x = order[i];
        while (j < n && (n - j) >= len[x]) add(sa.rank[j++]);
        res[x] = sa.sa[kth(k[x])] + 1;
    }

    for (i = 0; i < q; i++) write(res[i]);
    return 0;
}
