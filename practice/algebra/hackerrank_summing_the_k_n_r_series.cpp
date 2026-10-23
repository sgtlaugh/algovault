/***
 *
 * Lagrange Interpolation
 * Value of a polynomial at any x from its values at x = 0, 1, ..., n - 1, modulo a prime
 *
 * Complexity:
 *   - Lagrange(terms, mod): O(n + log mod)
 *   - interpolate: O(n) time and memory per call
 *   - get_terms(terms, mod, l, r): O(n (r - l + 1) + log mod)
 *   - find_degree: O(n (d + 2)) for the degree d it finds, O(n^2) when it returns -1
 *
 * Lagrange lagrange(terms, mod): the polynomial f of degree below n with f(x) = terms[x] for x = 0..n - 1
 * lagrange.interpolate(x): f(x) in [0, mod), for any long long x, negatives included
 * get_terms(terms, mod, l, r): f(l), f(l + 1), ..., f(r), for any ints l <= r
 * find_degree(terms, mod): the lowest degree d of a polynomial through all the terms, or -1 when fewer than d + 2 terms leave it unconfirmed
 *
 * mod must be a prime below 2^31 with n <= mod, which is asserted
 * Terms may be any ints, negatives included, and are reduced modulo mod
 * The zero polynomial has degree 0, so find_degree needs at least 2 terms of it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Lagrange{
    vector<int> terms, weight;
    int mod, n;

    Lagrange(const vector<int>& terms, int mod) : terms(terms), mod(mod){
        n = terms.size();
        assert(n <= mod);
        for (auto& t : this->terms){
            if ((t %= mod) < 0) t += mod;
        }

        int fact = 1;
        for (int i = 1; i < n; i++) fact = (long long)fact * i % mod;

        vector<int> inv_fact(n, expo(fact, mod - 2));
        for (int i = n - 1; i > 0; i--) inv_fact[i - 1] = (long long)inv_fact[i] * i % mod;

        /// The basis denominator, the product of (i - j) over j != i, is (-1)^(n - 1 - i) i! (n - 1 - i)!
        weight.resize(n);
        for (int i = 0; i < n; i++){
            weight[i] = (long long)inv_fact[i] * inv_fact[n - i - 1] % mod * this->terms[i] % mod;
        }
    }

    int interpolate(long long x){
        if (0 <= x && x < n) return terms[x];
        x %= mod;

        /// prefix[i] * suffix[i] is the product of (x - j) over j != i
        vector<int> prefix(n, 1), suffix(n, 1);
        for (int i = 1; i < n; i++){
            prefix[i] = (long long)prefix[i - 1] * (x - i + 1) % mod;
            if (prefix[i] < 0) prefix[i] += mod;
        }

        for (int i = n - 2; i >= 0; i--){
            suffix[i] = (long long)suffix[i + 1] * (x - i - 1) % mod;
            if (suffix[i] < 0) suffix[i] += mod;
        }

        long long res = 0;
        for (int i = 0; i < n; i++){
            int w = (long long)prefix[i] * suffix[i] % mod * weight[i] % mod;
            if ((n - 1 - i) & 1) w = mod - w;
            res += w;
        }

        return res % mod;
    }

    int expo(int a, int b){
        int res = 1;

        while (b){
            if (b & 1) res = (long long)res * a % mod;
            a = (long long)a * a % mod;
            b >>= 1;
        }

        return res;
    }
};

vector<int> get_terms(const vector<int>& terms, int mod, int l, int r){
    Lagrange lagrange(terms, mod);

    vector<int> res;
    for (int i = l; i <= r; i++){
        res.push_back(lagrange.interpolate(i));
    }

    return res;
}

int find_degree(const vector<int>& terms, int mod){
    int n = terms.size();
    assert(n <= mod);

    vector<int> diff(terms);
    for (auto& x : diff){
        if ((x %= mod) < 0) x += mod;
    }

    /// The k-th differences all vanish exactly when k > d, as the d-th ones equal d! times the leading coefficient and d < n <= mod
    /// Exact but slower than matching prefix fits at a few points, which roots at those points fool: 509 ms against 107 ms at n = 30000, d = n - 1
    for (int k = 0; k < n; k++){
        if (all_of(diff.begin(), diff.end() - k, [](int x){ return x == 0; })){
            int d = max(k - 1, 0);
            return n >= d + 2 ? d : -1;
        }

        for (int i = 0; i + k + 1 < n; i++){
            diff[i] = diff[i + 1] - diff[i];
            if (diff[i] < 0) diff[i] += mod;
        }
    }

    return -1;
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

#define MOD 1000000007

long long expo(long long a, long long b){
    long long res = 1;
    a %= MOD;
    while (b){
        if (b & 1) res = res * a % MOD;
        a = a * a % MOD;
        b >>= 1;
    }
    return res;
}

vector<int> powers(int d){
    int i, j;
    vector<int> pw(d + 2, 0), primes;
    vector<char> composite(d + 2, 0);

    pw[0] = (d == 0), pw[1] = 1;
    for (i = 2; i <= (d + 1); i++){
        if (!composite[i]) primes.push_back(i), pw[i] = expo(i, d);
        for (j = 0; j < (int)primes.size() && (long long)i * primes[j] <= (d + 1); j++){
            composite[i * primes[j]] = 1;
            pw[i * primes[j]] = (long long)pw[i] * pw[primes[j]] % MOD;
            if (i % primes[j] == 0) break;
        }
    }
    return pw;
}

long long solve(long long r, int d, long long n){
    int i;
    if (!r) return (n >= 1 && d == 0);

    vector<int> s = powers(d);
    long long cur = 1, sum = 0;
    for (i = 0; i <= (d + 1); i++){
        long long term = cur * s[i] % MOD;
        s[i] = sum;
        sum = (sum + term) % MOD;
        cur = cur * r % MOD;
    }

    if (r == 1){
        Lagrange lagrange(s, MOD);
        return lagrange.interpolate(n);
    }

    vector<int> inv(d + 3, 1);
    for (i = 2; i <= (d + 2); i++) inv[i] = (MOD - (long long)(MOD / i) * inv[MOD % i] % MOD) % MOD;

    long long ri = expo(r, MOD - 2), binom = 1, num = 0;
    for (i = 0, cur = 1; i <= (d + 1); i++){
        long long term = binom * s[i] % MOD * cur % MOD;
        if ((d + 1 - i) & 1) num = (num - term + MOD) % MOD;
        else num = (num + term) % MOD;
        binom = binom * (d + 1 - i) % MOD * inv[i + 1] % MOD;
        cur = cur * ri % MOD;
    }

    long long c = num * expo(expo((ri - 1 + MOD) % MOD, d + 1), MOD - 2) % MOD;

    s.pop_back();
    for (i = 0, cur = 1; i <= d; i++){
        s[i] = (s[i] - c + MOD) % MOD * cur % MOD;
        cur = cur * ri % MOD;
    }

    Lagrange lagrange(s, MOD);
    return (c + expo(r, n) * lagrange.interpolate(n)) % MOD;
}

int main(){
    int t, k;
    long long n, r;

    read(t);
    while (t--){
        read(k, n, r);
        write(solve(r % MOD, k, n + 1));
    }
    return 0;
}
