/***
 *
 * Fraction
 * Exact rational arithmetic in long long, always kept reduced with a positive denominator
 *
 * Complexity: O(log v) per arithmetic operation for the gcd, v the largest numerator or denominator, O(1) per comparison
 *
 * Fraction f(x): the integer x, 0 by default
 * Fraction f(num, den): num / den in lowest terms, den != 0
 * f.num, f.den: the reduced numerator and denominator, den > 0, so 1/-2 and -2/4 both hold -1 and 2
 * f + g, f - g, f * g: the sum, difference and product
 * f / g: the quotient, g != 0
 * -f and f += g, -=, *=, /=: negation and the compound assignments
 * f == g, !=, <, <=, >, >=: exact comparison
 *
 * Every binary operator and comparison also takes a long long on either side
 * Intermediate products use __int128, every reduced result must fit in long long (asserted)
 * Requires __int128 (64-bit GCC or Clang)
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

    /// Non-members so an integer works on either side, 1 + half as well as half + 1
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

        /// 64-bit gcd is about 2x faster, -LLONG_MAX keeps abs(LLONG_MIN) out of std::gcd
        __int128 g = (n >= -LLONG_MAX && n <= LLONG_MAX && d <= LLONG_MAX) ? __gcd(abs((long long)n), (long long)d) : gcd128(n, d);
        n /= g, d /= g;
        assert(n >= LLONG_MIN && n <= LLONG_MAX && d <= LLONG_MAX);

        Fraction f;
        f.num = (long long)n, f.den = (long long)d;
        return f;
    }

    static __int128 gcd128(__int128 a, __int128 b){
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        while (b){ __int128 r = a % b; a = b, b = r; }
        return a;
    }
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
 *
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ... and returns how many it read before the input ended
 * read accepts integers of any width including __int128, float, double, long double, strings, chars and vectors of these
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... separated by single spaces and ends the line, a vector argument on its own line
 * write accepts integers, strings, chars, vectors of these and fio::fixed values
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 * fio::flush(): writes the pending output now, which also happens at normal exit
 *
 * Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 * A char reads one non-space character
 * A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 * Writing a float without fixed() does not compile, so its precision is always chosen
 * Floats print at printf's speed, since the formatting itself dominates; a floating token keeps its first 127 characters
 * Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: stdin and stdout are one stream each per process, so a second reader instance
/// would split the input between two buffers and a second writer would reorder the output
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

    /// Integers: takes every type the plain overloads below miss, so __int128 needs no is_integral specialization, which breaks C++20
    template <typename T>
    bool read_one(T &x){
        int c = ' ', neg = 0;
        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

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

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
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

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
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

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

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
}

using namespace fio;

#define MAXE 105
#define MAXM 25

int rows, cols;
map<string, int> id;
Fraction mat[MAXE][MAXM];
string s;

int element(const string& m, int& p){
    string name(1, m[p++]);
    if (p < (int)m.size() && islower(m[p])) name += m[p++];
    if (!id.count(name)) id[name] = rows++;
    return id[name];
}

int number(const string& m, int& p){
    int x = 0;
    if (p >= (int)m.size() || !isdigit(m[p])) return 1;
    while (p < (int)m.size() && isdigit(m[p])) x = x * 10 + (m[p++] - '0');
    return x;
}

void parse(const string& m, int sign){
    int p = 0;
    while (p < (int)m.size()){
        if (m[p] != '('){
            int e = element(m, p);
            mat[e][cols] += sign * number(m, p);
            continue;
        }

        vector<pair<int, int>> inner;
        p++;
        while (m[p] != ')'){
            int e = element(m, p);
            inner.push_back(make_pair(e, number(m, p)));
        }
        p++;

        int k = number(m, p);
        for (auto& e : inner) mat[e.first][cols] += sign * e.second * k;
    }
    cols++;
}

void side(const string& t, int sign){
    int i, start = 0;
    for (i = 0; i <= (int)t.size(); i++){
        if (i == (int)t.size() || t[i] == '+'){
            parse(t.substr(start, i - start), sign);
            start = i + 1;
        }
    }
}

vector<long long> solve(){
    int i, j, k, rank = 0;
    vector<int> where(cols, -1);

    for (j = 0; j < cols && rank < rows; j++){
        for (i = rank; i < rows && mat[i][j] == 0; i++) {}
        if (i == rows) continue;
        for (k = 0; k < cols; k++) swap(mat[i][k], mat[rank][k]);

        Fraction pivot = mat[rank][j];
        for (k = 0; k < cols; k++) mat[rank][k] /= pivot;
        for (i = 0; i < rows; i++){
            if (i == rank || mat[i][j] == 0) continue;
            Fraction f = mat[i][j];
            for (k = 0; k < cols; k++) mat[i][k] -= f * mat[rank][k];
        }
        where[j] = rank++;
    }

    int free_col = 0;
    while (where[free_col] != -1) free_col++;

    vector<Fraction> x(cols, 0);
    x[free_col] = 1;
    for (j = 0; j < cols; j++){
        if (where[j] != -1) x[j] = -mat[where[j]][free_col];
    }

    long long lcm = 1, g = 0;
    for (j = 0; j < cols; j++) lcm = lcm / __gcd(lcm, x[j].den) * x[j].den;

    vector<long long> res(cols);
    for (j = 0; j < cols; j++){
        res[j] = abs(x[j].num) * (lcm / x[j].den);
        g = __gcd(g, res[j]);
    }
    for (j = 0; j < cols; j++) res[j] /= g;
    return res;
}

int main(){
    int i, j;

    while (read_line(s)){
        if (s.empty()) continue;
        if (s == ".") break;

        id.clear();
        rows = cols = 0;
        for (i = 0; i < MAXE; i++){
            for (j = 0; j < MAXM; j++) mat[i][j] = 0;
        }

        int arrow = s.find("->");
        side(s.substr(0, arrow), 1);
        side(s.substr(arrow + 2, s.size() - arrow - 3), -1);

        write(solve());
    }
    return 0;
}
