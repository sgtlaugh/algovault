/***
 *
 * Pell Equation
 * Integer solutions of x^2 - d y^2 = 1 and x^2 - d y^2 = -1, exact or modulo m
 *
 * Complexity: O(r + log k) operations on T, for r = O(sqrt(d) log d) the period of the continued fraction of sqrt(d)
 *
 * negative_pell<T>(d): smallest solution (x, y) of x^2 - d y^2 = -1 with y > 0, nullopt if there is none
 * pell<T>(d): smallest solution (x1, y1) of x^2 - d y^2 = 1 with y1 > 0
 * pell<T>(d, k): (x_k, y_k) with x_k + y_k sqrt(d) = (x1 + y1 sqrt(d))^k, for k >= 0
 * pell(d, k, m): (x_k mod m, y_k mod m) as long long, for k >= 0 and 1 <= m < 2^63
 *
 * Notes:
 *   - T is long long by default, any integer type with +, *, % and construction from long long works, like Bignum
 *   - Exact while the answer fits T: x1 fits long long for every d <= 276, __int128 for every d <= 1020
 *   - The (x_k, y_k) for k >= 1 are all the solutions of x^2 - d y^2 = 1 with x, y > 0, in increasing order
 *   - The solutions of x^2 - d y^2 = -1 are the odd powers of negative_pell(d), there are some iff r is odd
 *   - d a perfect square: pell gives (1, 0) for every k, negative_pell gives nullopt, except (0, 1) for d = 1
 *   - 1 <= d < 2^63, but r bounds the speed and x1 has O(r log d) digits
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Exact when m = 0
template<typename T>
T pell_reduce(const T& x, long long m){
    return m ? T(x % m) : x;
}

/// Last convergent p / q of the first period of sqrt(d), p^2 - d q^2 is -1 if the period is odd, else 1
template<typename T>
tuple<T, T, bool> pell_convergent(long long d, long long m){
    assert(d >= 1);
    long long a0 = sqrt((double)d);
    while (a0 * a0 > d) a0--;
    while ((__int128)(a0 + 1) * (a0 + 1) <= d) a0++;
    if (a0 * a0 == d) return make_tuple(pell_reduce(T(1), m), T(0), false);

    T p0 = pell_reduce(T(1), m), q0 = T(0), p = pell_reduce(T(a0), m), q = pell_reduce(T(1), m);
    long long num = 0, den = 1, a = a0;
    bool odd = false;
    while (true){
        num = den * a - num;
        den = (d - num * num) / den;
        a = (a0 + num) / den;
        odd = !odd;
        if (a == 2 * a0) break;

        T p1 = pell_reduce(a * p + p0, m), q1 = pell_reduce(a * q + q0, m);
        p0 = p, q0 = q, p = p1, q = q1;
    }

    return make_tuple(p, q, odd);
}

/// With m < 2^63 every value below is under 2^127, so unsigned __int128 never overflows
template<typename T>
pair<T, T> pell_multiply(const pair<T, T>& u, const pair<T, T>& v, long long d, long long m){
    T dyy = pell_reduce(pell_reduce(d * u.second, m) * v.second, m);
    return {pell_reduce(u.first * v.first + dyy, m), pell_reduce(u.first * v.second + u.second * v.first, m)};
}

template<typename T>
pair<T, T> pell_power(long long d, long long k, long long m){
    assert(k >= 0);
    T p, q;
    bool odd;
    tie(p, q, odd) = pell_convergent<T>(d, m);
    pair<T, T> base(p, q), res(pell_reduce(T(1), m), T(0));
    if (odd) base = pell_multiply(base, base, d, m);

    /// The last squaring is skipped, it would overflow an exact T even when the answer fits
    while (k){
        if (k & 1) res = pell_multiply(res, base, d, m);
        k >>= 1;
        if (k) base = pell_multiply(base, base, d, m);
    }

    return res;
}

template<typename T = long long>
pair<T, T> pell(long long d, long long k = 1){
    return pell_power<T>(d, k, 0);
}

pair<long long, long long> pell(long long d, long long k, long long m){
    assert(m >= 1);
    pair<unsigned __int128, unsigned __int128> xy = pell_power<unsigned __int128>(d, k, m);  /// unsigned measured 9% faster
    return {(long long)xy.first, (long long)xy.second};
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

const __int128 LIMIT = (__int128)16 * 1000000000000000000LL * 10000000000000000LL;

int squarefree(int c){
    for (int i = 3; (i * i) <= c; i += 2){
        while (c % (i * i) == 0) c /= (i * i);
    }

    return c;
}

int main(){
    int c, s, k;
    set<__int128> res;

    for (c = 1; c < 16; c += 2){
        s = squarefree(c);
        for (k = 1; ; k++){
            __int128 y = pell<__int128>(8 * s, k).second;
            __int128 x = (__int128)c * s * y * y;
            if (x >= LIMIT) break;
            res.insert(x);
        }
    }

    for (__int128 x : res) write(x);
    return 0;
}
