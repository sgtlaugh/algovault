/***
 *
 * Extended GCD
 * Extended gcd, modular inverse, and linear diophantine equations (solve and count solutions in a box)
 *
 * Complexity: O(log m) for extended_gcd, mod_inverse, diophantine and count_diophantine
 *
 * extended_gcd(a, b, x, y): returns g = gcd(a, b) and sets a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m for any 1 <= m < 2^63, -1 if gcd(a, m) != 1, a may be negative or >= m
 * diophantine(a, b, c, x, y): solves a * x + b * y = c, returns false if there is no solution
 *     a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
 * count_diophantine(a, b, c, x1, x2, y1, y2): number of solutions with x1 <= x <= x2 and y1 <= y <= y2
 *     every argument within [-1e9, 1e9], empty ranges give 0
 *     a = b = 0 and c = 0 gives the full box, (2e9 + 1)^2 at most, still within long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED extended_gcd
long long extended_gcd(long long a, long long b, long long& x, long long& y){
    assert(a >= 0 && b >= 0);

    long long x1 = 0, y1 = 1;
    x = 1, y = 0;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }

    return a;
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    long long x, y;
    if (extended_gcd(a, m, x, y) != 1) return -1;
    x %= m;
    return x < 0 ? x + m : x;
}
/// END SHARED extended_gcd

bool diophantine(long long a, long long b, long long c, long long& x, long long& y){
    if (a == 0 && b == 0){
        x = y = 0;
        return c == 0;
    }

    long long g = extended_gcd(llabs(a), llabs(b), x, y);
    if (c % g) return false;

    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;
}

long long count_diophantine(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    long long x, y;
    if (!diophantine(a, b, c, x, y)) return 0;
    if (a == 0 && b == 0) return max(0LL, x2 - x1 + 1) * max(0LL, y2 - y1 + 1);

    /// shrinks [lo, hi] to the k with lo_val <= base + k * step <= hi_val
    long long g = gcd(llabs(a), llabs(b)), lo = LLONG_MIN, hi = LLONG_MAX;
    auto restrict_steps = [&](long long base, long long step, long long lo_val, long long hi_val){
        if (step == 0){
            if (base < lo_val || base > hi_val) lo = 1, hi = 0;
            return;
        }

        if (step < 0) step = -step, base = -base, tie(lo_val, hi_val) = make_pair(-hi_val, -lo_val);
        auto floor_div = [&](long long num){ return num / step - (num % step < 0); };
        lo = max(lo, -floor_div(base - lo_val));
        hi = min(hi, floor_div(hi_val - base));
    };
    restrict_steps(x, b / g, x1, x2);
    restrict_steps(y, -a / g, y1, y2);
    return max(0LL, hi - lo + 1);
}

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

long long apply(long long x, char op, long long y, long long n){
    if (op == '+') return (x + y) % n;
    if (op == '-') return (x - y + n) % n;
    if (op == '*') return (unsigned __int128)x * y % n;

    long long inv = mod_inverse(y, n);
    if (inv == -1) return -1;
    return (unsigned __int128)x * inv % n;
}

int main(){
    long long n, x, y;
    int t;
    char op;

    while (read(n, t) == 2 && (n || t)){
        while (t--){
            read(x, op, y);
            write(apply(x, op, y, n));
        }
    }

    return 0;
}
