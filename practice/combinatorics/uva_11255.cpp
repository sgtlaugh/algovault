/***
 *
 * Necklaces
 * Counts rings of n beads in k colors up to rotation (necklaces), or up to rotation and turning over (bracelets)
 *
 * Complexity:
 *   - Any number of each color: O(sqrt(n) + d(n) log n) for d(n) divisors of n, at most 6720 for n <= 10^12
 *   - Fixed counts: O(sqrt(g) + d(g) (n + k)) for k colors, g the gcd of the counts
 *
 * necklaces<T>(n, k): number of necklaces, exact in T
 * bracelets<T>(n, k): number of bracelets, exact in T
 * necklaces(n, k, m): number of necklaces modulo m, needs n m < 2^64
 * bracelets(n, k, m): number of bracelets modulo m, needs 2 n m < 2^64
 * necklaces<T>(counts): number of necklaces with exactly counts[j] beads of color j, exact in T
 * bracelets<T>(counts): number of bracelets with exactly counts[j] beads of color j, exact in T
 *
 * Notes:
 *   - n >= 1, k >= 0 and m >= 1, and m may be any modulus, prime or not, even one dividing n
 *   - T is any integer type with +, * and / by a long long: long long, __int128, Bignum from bignum.cpp
 *   - T must hold n times the necklace count and 2 n times the bracelet count
 *   - counts: n beads in all, n >= 1, and T must hold 2 n times the n! / (counts[0]! counts[1]! ...) arrangements
 *   - counts modulo m: reduce the Bignum result, since the multinomials need exact division
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace burnside{
    /// A residue modulo mod = n m, whose division by n is exact on the integers and gives the answer modulo m
    struct Mod{
        unsigned long long v, mod;

        Mod operator + (const Mod& o) const{
            return {v >= mod - o.v ? v - (mod - o.v) : v + o.v, mod};
        }

        Mod operator * (const Mod& o) const{
            return {(unsigned long long)((unsigned __int128)v * o.v % mod), mod};
        }

        Mod operator * (long long x) const{
            return *this * Mod{(unsigned long long)x % mod, mod};
        }

        Mod operator / (long long x) const{
            return {v / x, mod};
        }
    };

    /// Left to right, so no intermediate exceeds b^e and an exact T never overflows before the answer does
    template<typename T>
    T power(const T& b, long long e){
        T res = b;
        for (int bit = __lg(e) - 1; bit >= 0; bit--){
            res = res * res;
            if (e >> bit & 1) res = res * b;
        }

        return res;
    }

    /// Every divisor d of n with phi(d)
    vector<pair<long long, long long>> divisors(long long n){
        vector<pair<long long, long long>> divs = {{1, 1}};
        long long rest = n;
        for (long long p = 2; rest > 1; p++){
            if (p > rest / p) p = rest;
            if (rest % p) continue;

            int len = divs.size();
            for (long long pj = 1; rest % p == 0; pj *= p){
                rest /= p;
                for (int i = 0; i < len; i++){
                    long long d = divs[i].first, phi = divs[i].second;
                    divs.push_back({d * pj * p, phi * pj * (p - 1)});
                }
            }
        }

        return divs;
    }

    /// (sum of c / len)! / prod (c / len)! over the counts c, rounding down, by exact divisions in T
    template<typename T>
    T multinomial(const vector<int>& counts, long long len){
        T res(1);
        long long total = 0;
        for (int c : counts){
            for (long long i = 1; i <= c / len; i++){
                total++;
                res = res * total / i;
            }
        }

        return res;
    }

    /// Sum of phi(d) k^(n / d) over the divisors d of n, n times the necklace count
    template<typename T>
    T rotations(long long n, const T& k){
        vector<pair<long long, long long>> divs = divisors(n);
        T res = power(k, n);
        for (size_t i = 1; i < divs.size(); i++) res = res + power(k, n / divs[i].first) * divs[i].second;
        return res;
    }

    /// Rotations with cycles of length d fix arrangements only when d divides every count, so d runs over g's divisors
    template<typename T>
    T rotations(const vector<int>& counts){
        long long g = 0;
        for (int c : counts){
            assert(c >= 0);
            g = __gcd(g, (long long)c);
        }
        assert(g >= 1);

        T res(0);
        for (auto x : divisors(g)) res = res + multinomial<T>(counts, x.first) * x.second;
        return res;
    }

    /// Sum of colorings fixed by each reflection, 2 n times the bracelet count minus rotations(n, k)
    template<typename T>
    T reflections(long long n, const T& k){
        if (n % 2) return power(k, (n + 1) / 2) * n;
        T half = power(k, n / 2);
        return (half + half * k) * (n / 2);
    }
}

template<typename T>
T necklaces(long long n, long long k){
    assert(n >= 1 && k >= 0);
    return burnside::rotations(n, T(k)) / n;
}

template<typename T>
T bracelets(long long n, long long k){
    assert(n >= 1 && k >= 0 && n < (1LL << 62));
    T kt(k);
    return (burnside::rotations(n, kt) + burnside::reflections(n, kt)) / (2 * n);
}

long long necklaces(long long n, long long k, long long m){
    assert(n >= 1 && k >= 0 && m >= 1 && (unsigned __int128)n * m <= ULLONG_MAX);
    unsigned long long mod = (unsigned long long)n * m;
    return (burnside::rotations(n, burnside::Mod{k % mod, mod}) / n).v;
}

long long bracelets(long long n, long long k, long long m){
    assert(n >= 1 && k >= 0 && m >= 1 && (unsigned __int128)2 * n * m <= ULLONG_MAX);
    unsigned long long mod = 2ULL * n * m;
    burnside::Mod kt{k % mod, mod};
    return ((burnside::rotations(n, kt) + burnside::reflections(n, kt)) / (2 * n)).v;
}

template<typename T>
T necklaces(const vector<int>& counts){
    T res = burnside::rotations<T>(counts);
    return res / accumulate(counts.begin(), counts.end(), 0LL);
}

/// Every reflection fixes n times multinomial(counts, 2) colorings in all, when at most 2 counts are odd, else none
template<typename T>
T bracelets(const vector<int>& counts){
    T res = burnside::rotations<T>(counts);
    long long n = accumulate(counts.begin(), counts.end(), 0LL), odd = 0;
    for (int c : counts) odd += c % 2;

    if (odd <= 2) res = res + burnside::multinomial<T>(counts, 2) * n;
    return res / (2 * n);
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
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - Reading into a char array does not compile: read a string, or a grid as vector<string>
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
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

    /// Without it a char * picks the integer template, an exact match that beats the const char * overload
    void write_one(char *s){
        write_one((const char *)s);
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
    int t, a, b, c;

    while (read(t)){
        while (t--){
            read(a, b, c);
            write((long long)bracelets<__int128>({a, b, c}));
        }
    }

    return 0;
}
