/***
 *
 * Matrix
 * Matrix products and powers modulo mod or over any semiring, plus transpose and rotation
 *
 * Complexity:
 *   - multiplication: O(n m p) for an n x m times an m x p matrix
 *   - pow(e): O(n^3 log e)
 *
 * Matrix a(n, m, mod): n x m zero matrix of long long modulo mod, 1 <= mod < 2^31
 * Matrix<T, Ops = PlusTimes<T>> a(n, m): n x m matrix of Ops::zero() with no modulus, like Matrix<long double> a(n, m)
 * Matrix<T, Ops>::identity(n, mod = 0): the n x n identity, written Matrix<>::identity(n, mod) when modular
 * a[i][j] = value: any T; with a modulus any long long, reduced modulo mod when multiplying
 * a * b: the product, needs a.m == b.n and the same mod
 * a.is_normalized(): whether every entry is already in [0, mod), always true without a modulus
 * a.normalized(): a copy with every entry in [0, mod), an unchanged copy without a modulus
 * a.pow(e): a^e for a square matrix and e >= 0, a.pow(0) is the identity
 * a.rotate(): the matrix turned 90 degrees clockwise
 * a.transpose(): the m x n transpose
 *
 * Notes:
 *   - Recurrence f(n) = 2 f(n - 1) + 3 f(n - 2): with M = {{2, 3}, {1, 0}}, (f(n + 1), f(n)) = M^n * (f(1), f(0))
 *   - Ops holds static zero(), one(), add(x, y), mul(x, y) of a semiring, so zero must absorb in mul
 *   - Min-plus: zero = INF, one = 0, add = min, mul = + capped at INF; a.pow(k)[u][v] is the shortest k-edge walk
 *   - Only the default Matrix<long long> takes a modulus, any other T or Ops asserts mod == 0
 *   - Without a modulus nothing is reduced: long long can overflow, unsigned long long wraps modulo 2^64
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct PlusTimes{
    static T zero(){ return 0; }
    static T one(){ return 1; }
    static T add(T x, T y){ return x + y; }
    static T mul(T x, T y){ return x * y; }
};

template<typename T = long long, typename Ops = PlusTimes<T>>
struct Matrix{
    static constexpr bool MODULAR = is_same<T, long long>::value && is_same<Ops, PlusTimes<T>>::value;

    int n, m;
    long long mod;
    vector<vector<T>> a;

    Matrix(int n, int m, long long mod = 0) : n(n), m(m), mod(mod), a(n, vector<T>(m, Ops::zero())){
        assert(mod == 0 || (MODULAR && 1 <= mod && mod < (1LL << 31)));
    }

    static Matrix identity(int n, long long mod = 0){
        Matrix res(n, n, mod);
        for (int i = 0; i < n; i++) res.a[i][i] = Ops::one();
        return res.normalized();
    }

    vector<T>& operator[](int i){
        return a[i];
    }

    const vector<T>& operator[](int i) const{
        return a[i];
    }

    Matrix operator*(const Matrix& other) const{
        assert(m == other.n && mod == other.mod);

        Matrix res(n, other.m, mod);
        if (other.is_normalized()) multiply(*this, other, res);
        else multiply(*this, other.normalized(), res);
        return res;
    }

    bool is_normalized() const{
        return true;
    }

    Matrix normalized() const{
        Matrix res = *this;

        return res;
    }

    Matrix pow(long long e) const{
        assert(n == m && e >= 0);

        /// tmp is the one scratch buffer, so the loop allocates nothing
        Matrix res = identity(n, mod), base = normalized(), tmp(n, n, mod);
        for (; e; e >>= 1){
            if (e & 1){
                multiply(res, base, tmp);
                swap(res, tmp);
            }
            multiply(base, base, tmp);
            swap(base, tmp);
        }

        return res;
    }

    Matrix rotate() const{
        Matrix res = transpose();
        for (auto& row : res.a) reverse(row.begin(), row.end());
        return res;
    }

    Matrix transpose() const{
        Matrix res(m, n, mod);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) res.a[j][i] = a[i][j];
        }

        return res;
    }

    /// res = x * y in res's existing rows; y must be normalized and res must not alias x or y
    static void multiply(const Matrix& x, const Matrix& y, Matrix& res){
        for (int i = 0; i < x.n; i++){
            vector<T>& acc = res.a[i];
            fill(acc.begin(), acc.end(), Ops::zero());
            for (int k = 0; k < x.m; k++){
                const T v = x.a[i][k];
                const vector<T>& row = y.a[k];
                for (int j = 0; j < y.m; j++) acc[j] = Ops::add(acc[j], Ops::mul(v, row[j]));
            }
        }
    }

    static void multiply_mod(const Matrix& x, const Matrix& y, Matrix& res){
        const long long mod = x.mod, limit = mod * mod;  /// sums stay below 2 * mod^2 < 2^63

        for (int i = 0; i < x.n; i++){
            vector<long long>& acc = res.a[i];
            fill(acc.begin(), acc.end(), 0);
            for (int k = 0; k < x.m; k++){
                long long v = x.a[i][k] % mod;
                if (v < 0) v += mod;
                if (!v) continue;
                const vector<long long>& row = y.a[k];
                for (int j = 0; j < y.m; j++){
                    acc[j] += v * row[j];
                    if (acc[j] >= limit) acc[j] -= limit;
                }
            }
            for (long long& value : acc) value %= mod;
        }
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
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
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

int main(){
    long long p, q, n;

    while (read(p, q, n) == 3){
        if (!n){
            write(2);
            continue;
        }

        Matrix<unsigned long long> step(2, 2), state(2, 1);
        step[0][0] = p, step[0][1] = -(unsigned long long)q, step[1][0] = 1;
        state[0][0] = p, state[1][0] = 2;
        write((long long)(step.pow(n - 1) * state)[0][0]);
    }

    return 0;
}
