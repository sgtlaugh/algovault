/***
 *
 * Ternary Search and Golden-Section Search
 * Optimum of a unimodal function: integer argmax, real argmin
 *
 * Complexity:
 *   - golden_section_min: O(log((hi - lo) / w)) calls of f for the final width w, about 1 per step
 *   - ternary_search_max: O(log(hi - lo)) calls of f, 2 per step
 *
 * golden_section_min(lo, hi, f, eps = 1e-9): x in [lo, hi] minimizing a real f, |lo|, |hi| <= 1e300
 * ternary_search_max(lo, hi, f): smallest i in [lo, hi] maximizing f(i), -4e18 <= lo <= hi <= 4e18, any comparable f
 *
 * Notes:
 *   - golden_section_min needs f strictly decreasing, then strictly increasing, with an optional flat bottom
 *   - golden_section_min stops once hi - lo <= eps * max(1, |lo| + |hi|), so keep eps >= 1e-12
 *   - Near a smooth minimum rounding limits x to about 1e-8 relative whatever eps is, f(x) is accurate to rounding
 *   - A kink at c is located to eps * max(1, |c|)
 *   - Maximum: golden_section_min on -f
 *   - ternary_search_max needs f(lo) < ... < f(p) >= f(p + 1) >= ... >= f(hi), plateaus only at or after the peak p
 *   - Plateaus before the peak: search g(i) = f(lo + hi - i), then i -> lo + hi - i is the largest argmax
 *   - Minimum: ternary_search_max on -f, which needs f strictly decreasing, then non-decreasing
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
double golden_section_min(double lo, double hi, F f, double eps = 1e-9){
    assert(lo <= hi);

    const double r = (sqrt(5.0) - 1) / 2;
    double x1 = hi - r * (hi - lo), x2 = lo + r * (hi - lo);
    double f1 = f(x1), f2 = f(x2);

    while (hi - lo > eps * max(1.0, abs(lo) + abs(hi))){
        /// Rounding drifts the reused point by 1/r per step until lo < x1 < x2 < hi breaks and the minimum is discarded
        if (f1 < f2){
            hi = x2, x2 = x1, f2 = f1;
            x1 = hi - r * (hi - lo), f1 = f(x1);
            double ideal = lo + r * (hi - lo);
            if (abs(x2 - ideal) > (hi - lo) / 64) x2 = ideal, f2 = f(x2);
        }
        else{
            lo = x1, x1 = x2, f1 = f2;
            x2 = lo + r * (hi - lo), f2 = f(x2);
            double ideal = hi - r * (hi - lo);
            if (abs(x1 - ideal) > (hi - lo) / 64) x1 = ideal, f1 = f(x1);
        }
    }

    return (lo + hi) / 2;
}

template<typename F>
long long ternary_search_max(long long lo, long long hi, F f){
    assert(lo <= hi);

    /// f(i) < f(i + 1) holds exactly for i < p, so this is a binary search for the first i where it fails
    while (lo < hi){
        long long mid = lo + (hi - lo) / 2;
        if (f(mid) < f(mid + 1)) lo = mid + 1;
        else hi = mid;
    }

    return lo;
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

const double PI = acos(-1.0);

int ax, ay, bx, by;

double angle(double x){
    double ux = ax - x, vx = bx - x;
    return atan2(abs(ux * by - vx * ay), ux * vx + (double)ay * by);
}

double neg_angle(double x){
    return -angle(x);
}

int main(){
    read(ax, ay, bx, by);
    if ((long long)ay * by <= 0){
        write(fixed(PI, 6));
        return 0;
    }

    double q = 0;
    if (ay != by) q = ax - (double)ay * (bx - ax) / (by - ay);

    double left = golden_section_min(q - 1e7, q, neg_angle, 1e-12);
    double right = golden_section_min(q, q + 1e7, neg_angle, 1e-12);
    write(fixed(max(angle(left), angle(right)), 6));
    return 0;
}
