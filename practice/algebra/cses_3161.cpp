/***
 *
 * GCD and LCM Convolution
 * c[k] = sum of a[i] * b[j] over gcd(i, j) = k, or over lcm(i, j) = k, exact or modulo mod
 *
 * Complexity:
 *   - DivisorLattice(n): O(n log log n)
 *   - each convolution or transform: O(m log log m) for vectors of size m + 1
 *
 * DivisorLattice lattice(n): serves any number of calls on vectors of size at most n + 1
 * lattice.gcd_convolution(a, b, mod = 0): c of size max(a.size(), b.size()) with c[k] = sum of a[i] * b[j] over gcd(i, j) = k
 * lattice.lcm_convolution(a, b, mod = 0): the same over lcm(i, j) = k, dropping pairs whose lcm is out of range
 * lattice.multiple_zeta(a, mod = 0): in place, a[d] becomes the sum of a[m] over the multiples m of d
 * lattice.multiple_mobius(a, mod = 0): in place, undoes multiple_zeta
 * lattice.divisor_zeta(a, mod = 0): in place, a[m] becomes the sum of a[d] over the divisors d of m
 * lattice.divisor_mobius(a, mod = 0): in place, undoes divisor_zeta
 *
 * Vectors are indexed by value from 1: the transforms ignore a[0] and the convolutions return c[0] = 0
 * mod = 0 is exact over long long, needing (sum |a[i]|) * (sum |b[j]|) < 2^63 for a convolution and sum |a[i]| < 2^63 for a transform
 * 1 <= mod < 2^62 gives results in [0, mod) for any long long input, negatives included, and mod need not be prime
 * Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DivisorLattice{
    int n;
    vector<int> primes;

    DivisorLattice(int n) : n(n){
        vector<bool> composite(n + 1, false);

        for (int i = 2; i <= n; i++){
            if (composite[i]) continue;
            primes.push_back(i);
            for (long long j = (long long)i * i; j <= n; j += i) composite[j] = true;
        }
    }

    vector<long long> gcd_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, true);
    }

    vector<long long> lcm_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, false);
    }

    void multiple_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) add(a[i], a[i * p], mod);
        }
    }

    void multiple_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) subtract(a[i], a[i * p], mod);
        }
    }

    void divisor_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) add(a[i * p], a[i], mod);
        }
    }

    void divisor_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) subtract(a[i * p], a[i], mod);
        }
    }

    vector<long long> convolve(vector<long long>& a, vector<long long>& b, long long mod, bool by_gcd) const{
        size_t size = max(a.size(), b.size());
        if (!size) return {};
        a.resize(size, 0), b.resize(size, 0);

        if (by_gcd) multiple_zeta(a, mod), multiple_zeta(b, mod);
        else divisor_zeta(a, mod), divisor_zeta(b, mod);

        for (size_t i = 1; i < size; i++){
            a[i] = mod ? (long long)((__int128)a[i] * b[i] % mod) : a[i] * b[i];
        }

        if (by_gcd) multiple_mobius(a, mod);
        else divisor_mobius(a, mod);
        a[0] = 0;
        return a;
    }

    int prepare(vector<long long>& a, long long mod) const{
        assert((long long)a.size() <= (long long)n + 1);
        if (mod) for (auto& x : a) x = (x % mod + mod) % mod;
        return (int)a.size() - 1;
    }

    static void add(long long& x, long long y, long long mod){
        x += y;
        if (mod && x >= mod) x -= mod;
    }

    static void subtract(long long& x, long long y, long long mod){
        x -= y;
        if (mod && x < 0) x += mod;
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

/// Because is_integral is not always true for __int128
#ifdef __SIZEOF_INT128__
    namespace std{
        template<> struct is_integral<__int128_t>: true_type {};
    }
#endif

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

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
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

    /// enable_if in the return type, since a second default template argument would redeclare the integral overload
    template <typename T>
    typename enable_if<is_floating_point<T>::value, bool>::type read_one(T& x){
        int c = ' ', len = 0;
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;

        /// Parsing at the target's own width, since strtold then a cast to double rounds twice and can miss by one ulp
        if (is_same<T, float>::value) x = strtof(tmpbuf, nullptr);
        else if (is_same<T, double>::value) x = strtod(tmpbuf, nullptr);
        else x = strtold(tmpbuf, nullptr);
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

    bool read_one(char& ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T>& v){
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

    int read(){
        return 0;
    }

    template <typename T, typename ...Args>
    int read(T &x, Args& ...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(const char& c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char* s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    void write_one(const string& s){
        for (auto &&c: s) write_char(c);
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
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

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

    void write_one(const Fixed& f){
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
    void write_one(const vector<T>& v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    void write() {}

    template <typename T> struct is_vector : false_type {};
    template <typename T, typename A> struct is_vector<vector<T, A>> : true_type {};

    template <typename T, typename ...Args>
    void write(const T& x, const Args& ...args){
        write_one(x);
        write_char(sizeof...(args) && !is_vector<T>::value ? ' ' : '\n');
        write(args...);
    }
}

using namespace fio;

#define MAX 200010
#define MOD 1000000007

long long pw[MAX];

int main(){
    int n, i, x;

    read(n);
    vector<long long> cnt(n + 1, 0);
    for (i = 0; i < n; i++){
        read(x);
        cnt[x]++;
    }

    pw[0] = 1;
    for (i = 1; i <= n; i++) pw[i] = pw[i - 1] * 2 % MOD;

    DivisorLattice lattice(n);
    lattice.multiple_zeta(cnt);
    for (i = 1; i <= n; i++) cnt[i] = pw[cnt[i]] - 1;
    lattice.multiple_mobius(cnt, MOD);

    cnt.erase(cnt.begin());
    write(cnt);
    return 0;
}
