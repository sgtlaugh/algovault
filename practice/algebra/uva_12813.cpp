/***
 *
 * Fraction
 * Exact rational arithmetic in long long, always kept reduced with a positive denominator
 *
 * Complexity: O(log v) per arithmetic operation, v the largest numerator or denominator, O(1) per comparison
 *
 * Fraction f(x): the integer x, 0 by default
 * Fraction f(num, den): num / den in lowest terms, den != 0
 * f + g, f - g, f * g: the sum, difference and product
 * f / g: the quotient, g != 0
 * -f and f += g, -=, *=, /=: negation and the compound assignments
 * f == g, !=, <, <=, >, >=: exact comparison
 * f.num, f.den: the reduced numerator and denominator, den > 0, so 1/-2 and -2/4 both hold -1 and 2
 *
 * Notes:
 *   - Every binary operator and comparison also takes a long long on either side
 *   - Every reduced result must fit in long long, intermediate products use __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Fraction{
    long long num, den;

    Fraction(long long x = 0) : num(x), den(1) {}

    Fraction(long long n, long long d){
        *this = make(n, d);
    }

    friend Fraction operator+(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.den + (__int128)b.num * a.den, (__int128)a.den * b.den);
    }

    friend Fraction operator-(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.den - (__int128)b.num * a.den, (__int128)a.den * b.den);
    }

    friend Fraction operator*(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.num, (__int128)a.den * b.den);
    }

    friend Fraction operator/(const Fraction& a, const Fraction& b){
        assert(b.num != 0);
        return make((__int128)a.num * b.den, (__int128)a.den * b.num);
    }

    friend Fraction operator-(const Fraction& a){
        return make(-(__int128)a.num, a.den);
    }

    Fraction& operator+=(const Fraction& b){
        return *this = *this + b;
    }

    Fraction& operator-=(const Fraction& b){
        return *this = *this - b;
    }

    Fraction& operator*=(const Fraction& b){
        return *this = *this * b;
    }

    Fraction& operator/=(const Fraction& b){
        return *this = *this / b;
    }

    friend bool operator==(const Fraction& a, const Fraction& b){
        return a.num == b.num && a.den == b.den;
    }

    friend bool operator!=(const Fraction& a, const Fraction& b){
        return !(a == b);
    }

    friend bool operator<(const Fraction& a, const Fraction& b){
        return (__int128)a.num * b.den < (__int128)b.num * a.den;
    }

    friend bool operator>(const Fraction& a, const Fraction& b){
        return b < a;
    }

    friend bool operator<=(const Fraction& a, const Fraction& b){
        return !(b < a);
    }

    friend bool operator>=(const Fraction& a, const Fraction& b){
        return !(a < b);
    }

    static Fraction make(__int128 n, __int128 d){
        assert(d != 0);
        if (d < 0) n = -n, d = -d;

        __int128 g = (n >= -LLONG_MAX && n <= LLONG_MAX && d <= LLONG_MAX) ? __gcd(llabs((long long)n), (long long)d) : gcd128(n, d);
        n /= g, d /= g;
        assert(n >= LLONG_MIN && n <= LLONG_MAX && d <= LLONG_MAX);

        Fraction f;
        f.num = (long long)n, f.den = (long long)d;
        return f;
    }

    static __int128 gcd128(__int128 a, __int128 b){
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        while (b) a %= b, swap(a, b);
        return a;
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

long long ten[20];

string solve(const string& str){
    int i = 0, n = str.size(), a = 0, b = 0;
    long long whole = 0, digits = 0, period = 0;
    bool neg = false;

    if (str[i] == '-' || str[i] == '+') neg = (str[i++] == '-');
    while (i < n && isdigit(str[i])) whole = whole * 10 + (str[i++] - '0');
    if (i < n && str[i] == '.') i++;
    while (i < n && isdigit(str[i])) digits = digits * 10 + (str[i++] - '0'), a++;
    if (i < n && str[i] == '('){
        i++;
        while (i < n && isdigit(str[i])) period = period * 10 + (str[i++] - '0'), b++;
    }

    Fraction part = Fraction(digits, ten[a]);
    if (b) part += Fraction(period, ten[a] * (ten[b] - 1));
    whole += part.num / part.den;
    long long rem = part.num % part.den;

    string sign = (neg && (whole || rem)) ? "-" : "";
    if (!rem) return sign + to_string(whole);
    if (!whole) return sign + to_string(rem) + "/" + to_string(part.den);
    return sign + to_string(whole) + " " + to_string(rem) + "/" + to_string(part.den);
}

int main(){
    int i;
    string str;

    for (ten[0] = 1, i = 1; i < 20; i++) ten[i] = ten[i - 1] * 10;

    while (read_line(str)){
        str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
        if (!str.empty()) write(solve(str));
    }

    return 0;
}
