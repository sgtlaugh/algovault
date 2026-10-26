/***
 *
 * Knuth Partition DP
 * Splits the prefix [0, i) into exactly k non-empty consecutive segments of minimum total cost
 *
 * Complexity: O(n^2 + n k) calls to cost whatever k is, O(n) memory
 *
 * knuth_partition_dp<T>(n, k, cost): dp[0..n] with dp[i] the cheapest split of the first i elements, k >= 0
 *
 * Notes:
 *   - cost(j, i) is the cost of the segment of elements j..i-1, called only for 0 <= j < i <= n
 *   - dp[i] = numeric_limits<T>::max() when i < k (fewer elements than segments)
 *   - Requires cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c) for a <= b <= c <= d
 *   - Such costs: (sum of a non-negative segment)^2, sum of non-negative pair weights inside it, convex in its length
 *   - Terms f(j) + g(i) of any sign may be added to any such cost
 *   - T is passed explicitly and fits the sum of k costs, like long long
 *   - For maximization, negate the cost and the answer
 *   - Same call, contract and result as divide_conquer_dp, so either can replace the other
 *   - Faster than divide_conquer_dp when k log2(n) > n, about 3x at k = n
 *   - Far slower for smaller k, since even k = 2 costs O(n^2) calls
 *   - A cost moving a window [j, i) by pointers makes about 2 moves per call, as in divide_conquer_dp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Cost>
vector<T> knuth_partition_dp(int n, int k, const Cost& cost){
    assert(k >= 0);

    const T inf = numeric_limits<T>::max();
    vector<T> dp(n + 1, inf);
    if (k > n) return dp;
    if (k == 0){
        dp[0] = 0;
        return dp;
    }

    for (int i = n; i >= 1; i--) dp[i] = cost(0, i);

    /// opt[n + 1] = n - 1 caps the first scan of every layer, layer 1 has opt[i] = 0
    vector<int> opt(n + 2, 0);
    opt[n + 1] = n - 1;
    for (int g = 2; g <= k; g++){
        /// i falls, so dp[j < i] and opt[i] still hold layer g - 1 while opt[i + 1] already holds layer g
        for (int i = n; i >= g; i--){
            int best = max(opt[i], g - 1), hi = min(opt[i + 1], i - 1);
            T best_val = dp[best] + cost(best, i);
            for (int j = best + 1; j <= hi; j++){
                T val = dp[j] + cost(j, i);
                if (val < best_val) best = j, best_val = val;
            }

            dp[i] = best_val;
            opt[i] = best;
        }
        dp[g - 1] = inf;
    }

    return dp;
}

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

#define MAX 1010

long long x[MAX], sum_w[MAX], sum_xw[MAX];

long long cost(int j, int i){
    return x[i - 1] * (sum_w[i] - sum_w[j]) - (sum_xw[i] - sum_xw[j]);
}

int main(){
    int n, k, i;
    long long w;

    while (read(n, k) == 2){
        for (i = 0; i < n; i++){
            read(x[i], w);
            sum_w[i + 1] = sum_w[i] + w;
            sum_xw[i + 1] = sum_xw[i] + x[i] * w;
        }

        write(knuth_partition_dp<long long>(n, k, cost)[n]);
    }

    return 0;
}
