/***
 *
 * Linear Sieve
 * Smallest prime factors, primes, Euler's phi, Mobius mu and any multiplicative function for 1..n in one pass
 *
 * Complexity: O(n) time, 13 bytes per entry plus 4 per prime, about 135 MB at n = 1e7
 *
 * LinearSieve sieve(n): tables for 0..n
 * sieve.primes: primes <= n in order
 * sieve.spf[x]: smallest prime factor of x >= 2
 * sieve.spf_exponent[x]: exponent of spf[x] in x >= 2
 * sieve.phi[x]: Euler's totient of x >= 1
 * sieve.mu[x]: Mobius function of x >= 1
 * sieve.factorize(x): prime factors of 1 <= x <= n with multiplicity, ascending, O(log x)
 * sieve.multiplicative(f): g[x] for 1 <= x <= n where g is multiplicative and g(p^k) = f(p, k, p^k), g[1] = 1, in O(n)
 *     e.g. divisor count: f = [](long long p, int k, long long pk){ return k + 1; }
 *
 * Notes:
 *   - n up to about 1e8 by memory, 1.3 GB; n = 2^31 - 1 overflows n + 1
 *   - multiplicative() returns a new table of 8 bytes per entry on every call
 *   - For mu alone, mobius_function.cpp needs 1 byte per entry and runs about 3x faster
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LinearSieve{
    int n;
    vector<int> spf, primes, phi, mu;
    vector<unsigned char> spf_exponent;  /// the exponent of spf[x] in x, 1 byte since it is at most 30

    LinearSieve(int n) : n(n), spf(n + 1, 0), phi(n + 1, 0), mu(n + 1, 0), spf_exponent(n + 1, 0){
        if (n >= 1) phi[1] = mu[1] = 1;

        for (int i = 2; i <= n; i++){
            if (!spf[i]){
                spf[i] = i, primes.push_back(i);
                phi[i] = i - 1, mu[i] = -1, spf_exponent[i] = 1;
            }
            for (int p : primes){
                long long x = (long long)i * p;
                if (p > spf[i] || x > n) break;
                spf[x] = p;
                if (p == spf[i]) phi[x] = phi[i] * p, mu[x] = 0, spf_exponent[x] = spf_exponent[i] + 1;
                else phi[x] = phi[i] * (p - 1), mu[x] = -mu[i], spf_exponent[x] = 1;
            }
        }
    }

    vector<int> factorize(int x) const{
        assert(1 <= x && x <= n);
        vector<int> res;
        for (; x > 1; x /= spf[x]) res.push_back(spf[x]);
        return res;
    }

    template <typename F>
    vector<long long> multiplicative(F f) const{
        vector<long long> g(n + 1, 0);

        if (n >= 1) g[1] = 1;
        for (int x = 2; x <= n; x++){
            /// Rebuilding p^k costs 15% against a stored power table (90 vs 77 ms at n = 1e7), which would cost 40 MB
            int p = spf[x], k = spf_exponent[x], pk = p;
            for (int i = 1; i < k; i++) pk *= p;
            g[x] = g[x / pk] * f((long long)p, k, (long long)pk);
        }

        return g;
    }
};

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

#define MAX 100010

using namespace fio;

int n, ar[MAX], cnt[MAX];
LinearSieve sieve(MAX);

vector<int> divisors(int x){
    int i, j, l, len, pk;
    vector<int> res = {1}, f = sieve.factorize(x);

    for (i = 0; i < (int)f.size(); i = j){
        len = res.size(), pk = 1;
        for (j = i; j < (int)f.size() && f[j] == f[i]; j++){
            pk *= f[i];
            for (l = 0; l < len; l++) res.push_back(res[l] * pk);
        }
    }

    return res;
}

int main(){
    int t, i;
    long long cur, res;

    read(t);
    while (t--){
        read(n);
        for (i = 0; i < n; i++) read(ar[i]);
        sort(ar, ar + n);

        res = 0;
        for (i = 0; i < n; i++){
            cur = 0;
            for (int d: divisors(ar[i])) cur += (long long)sieve.phi[d] * cnt[d]++;
            res += cur * (n - 1 - i);
        }

        memset(cnt, 0, sizeof(int) * (ar[n - 1] + 1));
        write(res);
    }
    return 0;
}
