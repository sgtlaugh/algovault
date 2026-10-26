/***
 *
 * Convex Hull Trick
 * Minimum or maximum of lines y = k * x + b at a point, for lines added in slope order
 *
 * Complexity: O(1) amortized per add_line and query_increasing, O(log n) per query
 *
 * ConvexHullTrick<T> hull: no lines yet, answers the minimum
 * ConvexHullTrick<T, true> hull: the same hull answering the maximum
 * hull.add_line(k, b): add y = k * x + b, k at most every earlier slope for the minimum, at least for the maximum
 * hull.query(x): the minimum over the lines at x, any x in any order, at least one line
 * hull.query_increasing(x): the same, for x never smaller than in the previous query_increasing call
 *
 * Notes:
 *   - Every k * x + b must fit in T, like |k|, |x| <= 1e9 and |b| <= 1e18 for long long, and needs __int128
 *   - T may be double or long double, every answer then within rounding of the exact optimum
 *   - Slopes in any order or lines only on part of the range: li_chao_tree.cpp
 *   - dp[i] = min over j < i of dp[j] + a[j] * x[i], a non-increasing: query(x[i]), then add_line(a[i], dp[i])
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, bool MAX = false>
struct ConvexHullTrick{
    /// The lower hull with strictly decreasing slopes; the maximum hull is stored negated
    vector<T> k, b;
    int ptr = 0;

    void add_line(T slope, T intercept){
        if (MAX) slope = -slope, intercept = -intercept;
        assert(k.empty() || slope <= k.back());
        if (!k.empty() && slope == k.back()){
            if (intercept >= b.back()) return;
            k.pop_back(), b.pop_back();
        }

        while (k.size() >= 2 && useless(k.size() - 2, k.size() - 1, slope, intercept)) k.pop_back(), b.pop_back();
        k.push_back(slope), b.push_back(intercept);
        ptr = min(ptr, (int)k.size() - 1);
    }

    T query(T x) const{
        assert(!k.empty());
        int lo = 0, hi = k.size() - 1;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (value(mid, x) >= value(mid + 1, x)) lo = mid + 1;
            else hi = mid;
        }

        return MAX ? -value(lo, x) : value(lo, x);
    }

    T query_increasing(T x){
        assert(!k.empty());
        while (ptr + 1 < (int)k.size() && value(ptr + 1, x) <= value(ptr, x)) ptr++;
        return MAX ? -value(ptr, x) : value(ptr, x);
    }

    T value(int i, T x) const{
        return k[i] * x + b[i];
    }

    /// Line j between lines i and new is never the minimum once new crosses i no later than j does
    bool useless(int i, int j, T slope, T intercept) const{
        /// A cast to __int128 would drop the fractions of floating-point lines
        if constexpr (is_floating_point_v<T>){
            return (long double)(intercept - b[i]) * (k[i] - k[j]) <= (long double)(b[j] - b[i]) * (k[i] - slope);
        }
        return (__int128)(intercept - b[i]) * (k[i] - k[j]) <= (__int128)(b[j] - b[i]) * (k[i] - slope);
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

#define MAX 200010

const double INF = 1e18;

long long sum[MAX];
double A[MAX], B[MAX];

int main(){
    int n, k, i, j, t;

    read(n, k);
    for (i = 0; i < n; i++){
        read(t);
        sum[i + 1] = sum[i] + t;
        A[i + 1] = A[i] + (double)sum[i + 1] / t;
        B[i + 1] = B[i] + 1.0 / t;
    }

    vector<double> dp(n + 1, INF), nxt(n + 1);
    dp[0] = 0;
    for (j = 1; j <= k; j++){
        ConvexHullTrick<double> hull;
        fill(nxt.begin(), nxt.begin() + j, INF);
        for (i = j; i <= n; i++){
            hull.add_line(-sum[i - 1], dp[i - 1] - A[i - 1] + sum[i - 1] * B[i - 1]);
            nxt[i] = A[i] + hull.query_increasing(B[i]);
        }
        swap(dp, nxt);
    }

    write(fixed(dp[n], 10));
    return 0;
}
