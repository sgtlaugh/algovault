/***
 *
 * Subsequence Automaton
 * Next occurrence of each symbol after any position, for greedy subsequence matching, and distinct subsequence counts
 *
 * Complexity:
 *   - SubsequenceAutomaton: O(n alphabet) build time and memory, O(1) per next, O(|t|) per is_subsequence
 *   - SparseSubsequenceAutomaton: O(n log n) build, O(log n) per next, O(|t| log n) per is_subsequence
 *   - distinct_subsequences: O(n log n)
 *
 * SubsequenceAutomaton sa(s, alphabet = 26, first = 'a'): string, C string or vector s over [first, first + alphabet)
 * sa.is_subsequence(t): true if t is a subsequence of s, the empty t included
 * sa.next(i, c): smallest j >= i with s[j] = c, n if none or c is outside the alphabet, for 0 <= i <= n + 1
 *
 * SparseSubsequenceAutomaton sa(s): vector s of any values with < and ==, such as ints up to 1e9
 * sa.is_subsequence(t): true if the vector t is a subsequence of s, the empty t included
 * sa.next(i, c): smallest j >= i with s[j] = c, n if none, for 0 <= i <= n + 1
 *
 * distinct_subsequences(s, mod = 998244353): distinct non-empty subsequences of s modulo mod, 1 <= mod <= 2^62
 *
 * Notes:
 *   - next(n + 1, c) = n, so a greedy walk i = next(i, c) + 1 never leaves the table and matched t iff i <= n
 *   - SubsequenceAutomaton takes 4 (n + 2) alphabet bytes, 104 MB at n = 1e6 and alphabet 26, else use the sparse one
 *   - Fit [first, first + alphabet) to the symbols s uses: t may hold others, and time and memory grow with alphabet
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SubsequenceAutomaton{
    int n, alphabet, first;
    vector<int> next_pos;

    template<typename Container>
    SubsequenceAutomaton(const Container& s, int alphabet = 26, int first = 'a')
        : n(s.size()), alphabet(alphabet), first(first), next_pos((n + 2) * alphabet, n){
        for (int i = n - 1; i >= 0; i--){
            long long c = s[i] - (long long)first;
            assert(0 <= c && c < alphabet);
            copy_n(next_pos.begin() + (i + 1) * alphabet, alphabet, next_pos.begin() + i * alphabet);
            next_pos[i * alphabet + c] = i;
        }
    }

    SubsequenceAutomaton(const char* s, int alphabet = 26, int first = 'a')
        : SubsequenceAutomaton(string(s), alphabet, first) {}

    template<typename Container>
    bool is_subsequence(const Container& t) const{
        int i = 0;

        for (long long c : t){
            if (i > n) return false;
            i = next(i, c) + 1;
        }

        return i <= n;
    }

    bool is_subsequence(const char* t) const{
        return is_subsequence(string(t));
    }

    /// long long so a large value is never narrowed into the alphabet
    int next(int i, long long c) const{
        c -= first;
        if (c < 0 || c >= alphabet) return n;

        return next_pos[i * alphabet + c];
    }
};

template<typename T>
struct SparseSubsequenceAutomaton{
    int n;
    vector<pair<T, int>> occurrences;

    SparseSubsequenceAutomaton(const vector<T>& s) : n(s.size()){
        for (int i = 0; i < n; i++) occurrences.emplace_back(s[i], i);
        sort(occurrences.begin(), occurrences.end());
    }

    bool is_subsequence(const vector<T>& t) const{
        int i = 0;

        for (const T& c : t){
            if (i > n) return false;
            i = next(i, c) + 1;
        }

        return i <= n;
    }

    int next(int i, const T& c) const{
        auto it = lower_bound(occurrences.begin(), occurrences.end(), make_pair(c, i));

        return it != occurrences.end() && it->first == c ? it->second : n;
    }
};

template<typename Container>
long long distinct_subsequences(const Container& s, long long mod = 998244353){
    assert(1 <= mod && mod <= (1LL << 62));
    int n = s.size();
    vector<int> order(n), prev(n, -1);

    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){ return make_pair(s[a], a) < make_pair(s[b], b); });
    for (int k = 1; k < n; k++){
        if (s[order[k]] == s[order[k - 1]]) prev[order[k]] = order[k - 1];
    }

    /// count[i] = distinct subsequences of s[0..i), the empty one included; those ending at s[i] were new after prev
    vector<long long> count(n + 1);
    count[0] = 1;
    for (int i = 0; i < n; i++){
        long long repeated = prev[i] >= 0 ? count[prev[i]] : 0;
        count[i + 1] = ((2 * count[i] - repeated) % mod + mod) % mod;
    }

    return (count[n] - 1 + mod) % mod;
}

long long distinct_subsequences(const char* s, long long mod = 998244353){
    return distinct_subsequences(string(s), mod);
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

#define MAX 105

const long long INF = 1e12;

long long cnt[MAX][MAX];

int main(){
    int n, i, l;
    long long k, res = 0;
    string s;

    read(n, k, s);
    SubsequenceAutomaton sa(s);

    for (i = n + 1; i >= 0; i--){
        cnt[i][0] = 1;
        if (i > n) continue;

        for (char c = 'a'; c <= 'z'; c++){
            int j = sa.next(i, c);
            if (j == n) continue;

            for (l = 1; l <= n; l++) cnt[i][l] = min(INF, cnt[i][l] + cnt[j + 1][l - 1]);
        }
    }

    for (l = n; l >= 0 && k; l--){
        long long x = min(k, cnt[0][l]);
        res += (x * (n - l)), k -= x;
    }

    write(k ? -1 : res);
    return 0;
}
