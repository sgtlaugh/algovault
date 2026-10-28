/***
 *
 * Necklaces
 * Counts rings of n beads in k colors up to rotation (necklaces), or up to rotation and turning over (bracelets)
 *
 * Complexity: O(sqrt(n) + d(n) log n) for d(n) divisors of n, at most 6720 for n <= 10^12
 *
 * necklaces<T>(n, k): number of necklaces, exact in T
 * bracelets<T>(n, k): number of bracelets, exact in T
 * necklaces(n, k, m): number of necklaces modulo m, needs n m < 2^64
 * bracelets(n, k, m): number of bracelets modulo m, needs 2 n m < 2^64
 *
 * Notes:
 *   - n >= 1, k >= 0 and m >= 1, and m may be any modulus, prime or not, even one dividing n
 *   - T is any integer type with +, * and / by a long long: long long, __int128, Bignum from bignum.cpp
 *   - T must hold n times the necklace count and 2 n times the bracelet count
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

    /// Sum of phi(d) k^(n / d) over the divisors d of n, n times the necklace count
    template<typename T>
    T rotations(long long n, const T& k){
        vector<pair<long long, long long>> divisors = {{1, 1}};
        long long rest = n;
        for (long long p = 2; rest > 1; p++){
            if (p > rest / p) p = rest;
            if (rest % p) continue;

            int len = divisors.size();
            for (long long pj = 1; rest % p == 0; pj *= p){
                rest /= p;
                for (int i = 0; i < len; i++){
                    long long d = divisors[i].first, phi = divisors[i].second;
                    divisors.push_back({d * pj * p, phi * pj * (p - 1)});
                }
            }
        }

        T res = power(k, n);
        for (size_t i = 1; i < divisors.size(); i++) res = res + power(k, n / divisors[i].first) * divisors[i].second;
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

const int MOD = 1000000007;

int main(){
    int t, n, m, T = 0;
    char buf[64];

    while (read(t)){
        while (t--){
            read(n, m);
            snprintf(buf, sizeof(buf), "Case #%d: %lld", ++T, bracelets(m, n, MOD));
            write(string(buf));
        }
    }

    return 0;
}
