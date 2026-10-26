/***
 *
 * Slope Trick
 * Maintains a convex piecewise linear function f(x) with integer slopes as two heaps of breakpoints
 *
 * Complexity:
 *   - add_*: O(log n) for n breakpoints added so far
 *   - shift, prefix_min, suffix_min, argmin, get_min: O(1)
 *   - eval: O(n)
 *
 * SlopeTrick<T> f: starts as f(x) = 0, T is a signed integer type
 * add_a_minus_x(a): f(x) += max(0, a - x)
 * add_abs(a): f(x) += |x - a|
 * add_const(c): f(x) += c
 * add_x_minus_a(a): f(x) += max(0, x - a)
 * prefix_min(): f(x) = min of f(y) over y <= x
 * shift(lo, hi): f(x) = min of f(y) over x - hi <= y <= x - lo, needs lo <= hi
 * shift(d): f(x) = f(x - d), moves the graph right by d
 * suffix_min(): f(x) = min of f(y) over y >= x
 * argmin(): [l, r], the range where f is minimal, -INF or INF when it is unbounded on that side, INF = max of T / 2
 * eval(x): f(x), without changing f
 * get_min(): the minimum of f
 *
 * Notes:
 *   - Breakpoints, each side's total shift and eval arguments stay in [-1e18, 1e18] for long long, [-1e9, 1e9] for int
 *   - The minimum and any evaluated value fit in T
 *   - prefix_min() resets the right side's total shift to 0, suffix_min() the left side's
 *   - min_cost_non_decreasing(a) makes a non-decreasing at minimum sum |a[i] - b[i]|, run it on a[i] - i for increasing
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct SlopeTrick{
    /// above every allowed breakpoint, so an empty heap's sentinel adds 0 to min_f, and INF - (-INF) still fits in T
    static constexpr T INF = numeric_limits<T>::max() / 2;

    T min_f = 0, add_left = 0, add_right = 0;
    vector<T> left, right;  /// max-heap and min-heap, each left breakpoint <= each right one, stored minus the lazy add

    void add_a_minus_x(T a){
        min_f += max(T(0), a - top_right());
        push_right(a);
        push_left(pop_right());
    }

    void add_abs(T a){
        add_x_minus_a(a);
        add_a_minus_x(a);
    }

    void add_const(T c){
        min_f += c;
    }

    void add_x_minus_a(T a){
        min_f += max(T(0), top_left() - a);
        push_left(a);
        push_right(pop_left());
    }

    void prefix_min(){
        right.clear();
        add_right = 0;
    }

    void shift(T lo, T hi){
        assert(lo <= hi);
        add_left += lo, add_right += hi;
    }

    void shift(T d){
        shift(d, d);
    }

    void suffix_min(){
        left.clear();
        add_left = 0;
    }

    pair<T, T> argmin() const{
        return {top_left(), top_right()};
    }

    T eval(T x) const{
        T res = min_f;
        for (T l : left) res += max(T(0), l + add_left - x);
        for (T r : right) res += max(T(0), x - r - add_right);
        return res;
    }

    T get_min() const{
        return min_f;
    }

    T top_right() const{
        return right.empty() ? INF : right[0] + add_right;
    }

    void push_right(T x){
        right.push_back(x - add_right);
        push_heap(right.begin(), right.end(), greater<T>());
    }

    void push_left(T x){
        left.push_back(x - add_left);
        push_heap(left.begin(), left.end());
    }

    T pop_right(){
        T x = top_right();
        pop_heap(right.begin(), right.end(), greater<T>());
        right.pop_back();
        return x;
    }

    T top_left() const{
        return left.empty() ? -INF : left[0] + add_left;
    }

    T pop_left(){
        T x = top_left();
        pop_heap(left.begin(), left.end());
        left.pop_back();
        return x;
    }
};

/// min of sum |a[i] - b[i]| over non-decreasing b
long long min_cost_non_decreasing(const vector<long long>& a){
    SlopeTrick<long long> f;
    for (long long x : a){
        f.prefix_min();
        f.add_abs(x);
    }
    return f.get_min();
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

int main(){
    int n, i;

    read(n);
    vector<long long> ar(n);
    read(ar);
    for (i = 0; i < n; i++) ar[i] -= i;

    write(min_cost_non_decreasing(ar));
    return 0;
}
