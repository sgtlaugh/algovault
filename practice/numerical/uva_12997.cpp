/***
 *
 * Numerical Integration
 * Definite integral of f over [a, b] with composite Simpson, adaptive Simpson or Romberg
 *
 * Complexity:
 *   - adaptive_simpson: O(panels it needs) evaluations, recursion depth <= 44
 *   - romberg: 2^k + 1 evaluations for k levels, 8 <= k <= max_levels
 *   - simpson: 2n + 1 evaluations
 *
 * adaptive_simpson(a, b, f, eps = 1e-9): splits each of 64 equal start panels until its halves agree with it
 * romberg(a, b, f, eps = 1e-9, max_levels = 20): Richardson extrapolation of trapezoid sums on 2^k subintervals
 * simpson(a, b, f, n = 1000): composite Simpson on 2n equal subintervals for n >= 1, exact for cubics
 *
 * Notes:
 *   - f must be finite at every point of [a, b], endpoints included, and a > b gives the negated integral
 *   - eps is absolute and judged by the method's own error estimate, so a coincidental cancellation stops early
 *   - On random integrands the error exceeded eps + 1e-13 relative in about 1 in 10^5 calls, by at most 6x
 *   - adaptive_simpson and romberg also stop at 1e-14 of the integral of |f|, so an eps below rounding still terminates
 *   - adaptive_simpson refines only around kinks, jumps and infinite-slope endpoints like sqrt at 0
 *   - adaptive_simpson can miss a spike, or two jumps, within one start panel of width (b - a) / 64
 *   - romberg reads periodic f as constant when (b - a) / 256 is a whole number of periods
 *   - romberg returns its last estimate after max_levels, it is fast only on smooth f: use adaptive_simpson for kinks
 *   - simpson errs by (b - a) h^4 max|f''''| / 180 with h = (b - a) / 2n
 *   - Nest calls for multiple integrals:
 *         adaptive_simpson(0, 1, [](double x){ return adaptive_simpson(0, x, [&](double y){ return x * y; }); });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
double adaptive_simpson_panel(F& f, double a, double b, double fa, double fm, double fb, double whole, double eps, double& noise, int depth){
    const int MAX_DEPTH = 44;
    double m = (a + b) / 2, flm = f((a + m) / 2), frm = f((m + b) / 2);
    double left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    double delta = left + right - whole;

    /// The start grid can miss where |f| is large, so each panel's integral of |f| lifts the floor over its rounding
    noise = max(noise, 1e-14 * abs(b - a) / 12 * (abs(fa) + 4 * abs(flm) + 2 * abs(fm) + 4 * abs(frm) + abs(fb)));
    if (depth >= MAX_DEPTH || abs(delta) <= max(15 * eps, noise)) return left + right + delta / 15;

    return adaptive_simpson_panel(f, a, m, fa, flm, fm, left, eps / 2, noise, depth + 1)
         + adaptive_simpson_panel(f, m, b, fm, frm, fb, right, eps / 2, noise, depth + 1);
}

template<typename F>
double adaptive_simpson(double a, double b, F f, double eps = 1e-9){
    const int PANELS = 64;
    double h = (b - a) / (2 * PANELS), abs_integral = 0;
    vector<double> y(2 * PANELS + 1);

    for (int i = 0; i <= 2 * PANELS; i++){
        y[i] = f(i == 2 * PANELS ? b : a + i * h);
        abs_integral += abs(h / 3 * y[i]) * (i == 0 || i == 2 * PANELS ? 1 : i % 2 ? 4 : 2);
    }

    /// Not halved per level like eps: a panel's rounding noise shrinks with its width
    double noise = 1e-14 * abs_integral, res = 0;
    for (int i = 0; i < 2 * PANELS; i += 2){
        double whole = h / 3 * (y[i] + 4 * y[i + 1] + y[i + 2]), r = i + 2 == 2 * PANELS ? b : a + (i + 2) * h;
        res += adaptive_simpson_panel(f, a + i * h, r, y[i], y[i + 1], y[i + 2], whole, eps / PANELS, noise, 0);
    }

    return res;
}

template<typename F>
double romberg(double a, double b, F f, double eps = 1e-9, int max_levels = 20){
    const int MIN_LEVELS = 8;
    double h = b - a, fa = f(a), fb = f(b), abs_trapezoid = abs(h / 2) * (abs(fa) + abs(fb));
    vector<double> prev = {h / 2 * (fa + fb)}, cur;

    for (int k = 1; k <= max_levels; k++){
        h /= 2;
        double sum = 0, abs_sum = 0;
        for (long long i = 1; i < (1LL << k); i += 2){
            double y = f(a + i * h);
            sum += y;
            abs_sum += abs(y);
        }

        cur.assign(k + 1, 0);
        cur[0] = prev[0] / 2 + h * sum;
        abs_trapezoid = abs_trapezoid / 2 + abs(h) * abs_sum;
        double power = 4;
        for (int j = 1; j <= k; j++, power *= 4) cur[j] = cur[j - 1] + (cur[j - 1] - prev[j - 1]) / (power - 1);

        if (k >= MIN_LEVELS && abs(cur[k] - prev[k - 1]) <= max(eps, 1e-14 * abs_trapezoid)) return cur[k];
        swap(prev, cur);
    }

    return prev.back();
}

template<typename F>
double simpson(double a, double b, F f, int n = 1000){
    double h = (b - a) / (2 * n), sum = f(a) + f(b);
    for (int i = 1; i < 2 * n; i++) sum += f(a + i * h) * (i & 1 ? 4 : 2);
    return sum * h / 3;
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

int x, d;

double walk(double s){
    double a = hypot(x, s), b = hypot(x, s + d), m = (d + b - a) / 2;
    return (m * m / 2 + a * m + (d - m) * (d - m) / 2 + b * (d - m)) / d;
}

int main(){
    int t, y1, y2, T;
    double res;
    char buf[64];

    while (read(t)){
        T = 0;
        while (t--){
            read(x, y1, y2, d);

            if ((y1 - y2) == d) res = walk(y2);
            else res = adaptive_simpson(y2, y1 - d, walk) / (y1 - d - y2);

            snprintf(buf, sizeof(buf), "Case %d: %.6f", ++T, res);
            write(string(buf));
        }
    }

    return 0;
}
