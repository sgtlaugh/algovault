/***
 *
 * Fast Fourier Transform
 * Polynomial multiplication, exact over integers or modulo any mod, and binary string matching
 *
 * Complexity:
 *   - every call: O((n + m) log(n + m)) for inputs of sizes n and m
 *   - the first call at a new transform length len, the smallest power of two >= n + m: O(len) more to extend the roots of unity
 *   - memory: up to 4 complex arrays of len entries per precision, kept between calls, 32 bytes each in long double, so ~270 MB at len = 2^21
 *
 * FFT fft: an instance that keeps its roots of unity and work buffers between calls
 * fft.multiply(a, b): the linear product of a and b, negatives allowed, exact whenever every coefficient fits in a long long
 * fft.square(a): the same as multiply(a, a)
 * fft.mod_multiply(a, b, mod): the product modulo mod, 1 <= mod < 2^31, any input values, results in [0, mod)
 * fft.ll_multiply(a, b): the exact product of values in [0, LL_MULTIPLY_LIMIT), every coefficient below 2^63
 * fft.hamming_distance(str, pat): mismatches of binary pat against each window of str, |str| - |pat| + 1 values
 * fft.common_ones(str, pat): common 1 bits of binary pat and each window of str, |str| - |pat| + 1 values
 *
 * Products have length n + m - 1, and an empty input gives an empty result
 * Once max(n, m) * max|value|^2 passes EXACT_LIMIT (~9.2e16), multiply and square need |values| < LL_MULTIPLY_LIMIT and take about 2x longer
 * mod_multiply takes that slower path unless max(n, m) * (mod - 1)^2 <= EXACT_LIMIT, which holds for mod below about 3e5 at n = 10^6
 * The string methods return an empty vector when |pat| > |str| and |str| + 1 zeros for an empty pat
 * Use ntt.cpp modulo an NTT-friendly prime like 998244353 (about 4x faster) or for values near 10^9 (about 2x), fft for small values
 *
 * Circular convolution of two size n vectors: fold the linear product, add % mod for mod_multiply
 *   auto c = fft.multiply(a, b);
 *   for (int i = n; i < (int)c.size(); i++) c[i - n] += c[i];
 *   c.resize(n);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FFT{
    template <typename T>
    struct Complex{
        T real, img;

        Complex() {}  /// left uninitialised so growing a work buffer does not write it twice, every entry is written before it is read
        Complex(T real, T img=0) : real(real), img(img) {}

        Complex conjugate() const{
            return Complex(real, -img);
        }

        Complex operator + (const Complex& x) const{
            return Complex(real + x.real, img + x.img);
        }

        Complex operator - (const Complex& x) const{
            return Complex(real - x.real, img - x.img);
        }

        Complex operator * (const Complex& x) const{
            return Complex(real * x.real - img * x.img, real * x.img + img * x.real);
        }
    };

    /// Roots of unity, bit reversal and work buffers of one precision, kept between calls
    template <typename T>
    struct Engine{
        int last_len = -1;
        vector<int> rev;

        /// Work buffers are reused across calls: allocating them per call measured 10-25% slower at len = 2^20
        vector<Complex<T>> u, w, f, g;

        /// roots[k + j] = e^(i * pi * j / k) for every power of two k < roots.size(), grown on demand
        vector<Complex<T>> roots = {Complex<T>(0), Complex<T>(1)};

        void prepare(int len, int roots_len){
            for (int k = roots.size(); k < roots_len; k <<= 1){
                T theta = acosl(-1.0L) / k;
                Complex<T> mul(cos(theta), sin(theta));

                roots.resize(2 * k);
                for (int j = k >> 1; j < k; j++){
                    roots[2 * j] = roots[j];
                    roots[2 * j + 1] = roots[j] * mul;
                }
            }

            if ((int)u.size() < len) u.resize(len), w.resize(len), f.resize(len), g.resize(len), rev.resize(len);
            if (len != last_len){
                last_len = len;
                int bits = __lg(len);
                for (int i = 1; i < len; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bits - 1));
            }
        }

        void transform(const vector<Complex<T>>& in, vector<Complex<T>>& out, int len){
            for (int i = 0; i < len; i++) out[i] = in[rev[i]];

            for (int k = 1; k < len; k <<= 1){
                for (int i = 0; i < len; i += 2 * k){
                    for (int j = 0; j < k; j++){
                        Complex<T> z = out[i + j + k] * roots[k + j];
                        out[i + j + k] = out[i + j] - z;
                        out[i + j] = out[i + j] + z;
                    }
                }
            }
        }
    };

    static constexpr int SPLIT_BITS = 15;
    static constexpr int LL_MULTIPLY_LIMIT = 1500000000;

    /// Largest n_max * max|value|^2 one double or long double transform keeps exact, at least 4x below the first error up to n = 10^6
    static constexpr double DOUBLE_LIMIT = 0.01 / numeric_limits<double>::epsilon();
    static constexpr long double EXACT_LIMIT = 0.01L / numeric_limits<long double>::epsilon();

    Engine<double> fast;
    Engine<long double> wide;

    vector<long long> multiply(vector<long long> v1, vector<long long> v2){
        if (v1.empty() || v2.empty()) return {};
        if (v1 == v2) return square(v1);
        return product(v1, v2, false);
    }

    vector<long long> square(vector<long long> v){
        if (v.empty()) return {};
        vector<long long> copy = v;
        return product(v, copy, true);
    }

    vector<long long> mod_multiply(vector<long long> v1, vector<long long> v2, int mod){
        if (v1.empty() || v2.empty()) return {};
        for (auto& x: v1) x = (x % mod + mod) % mod;
        for (auto& x: v2) x = (x % mod + mod) % mod;

        /// Reduced values below a small mod keep the whole product exact in one transform, reduced afterwards
        if ((long double)max(v1.size(), v2.size()) * (mod - 1) * (mod - 1) <= EXACT_LIMIT){
            vector<long long> res = product(v1, v2, v1 == v2);
            for (auto& x: res) x %= mod;
            return res;
        }

        if (split_fits_double(v1, v2)) return split_mod(fast, v1, v2, mod);
        return split_mod(wide, v1, v2, mod);
    }

    vector<long long> ll_multiply(vector<long long> v1, vector<long long> v2){
        for (auto x: v1) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        for (auto x: v2) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        if (v1.empty() || v2.empty()) return {};
        return product(v1, v2, v1 == v2);
    }

    /// Mapping bits to +1 and -1 makes each window's dot product matches - mismatches
    vector<long long> hamming_distance(const string& str, const string& pat){
        assert(is_binary_string(str) && is_binary_string(pat));
        int n = str.size(), m = pat.size();
        if (m == 0) return vector<long long>(n + 1, 0);

        vector<long long> v1(n), v2(m);
        for (int i = 0; i < n; i++) v1[i] = str[i] == '1' ? 1 : -1;
        for (int i = 0; i < m; i++) v2[i] = pat[m - 1 - i] == '1' ? 1 : -1;

        auto dot = multiply(v1, v2);
        vector<long long> res;
        for (int i = 0; i + m <= n; i++) res.push_back((m - dot[i + m - 1]) / 2);
        return res;
    }

    vector<long long> common_ones(const string& str, const string& pat){
        assert(is_binary_string(str) && is_binary_string(pat));
        int n = str.size(), m = pat.size();
        if (m == 0) return vector<long long>(n + 1, 0);

        vector<long long> v1(n), v2(m);
        for (int i = 0; i < n; i++) v1[i] = str[i] - '0';
        for (int i = 0; i < m; i++) v2[i] = pat[m - 1 - i] - '0';

        auto dot = multiply(v1, v2);
        vector<long long> res;
        for (int i = 0; i + m <= n; i++) res.push_back(dot[i + m - 1]);
        return res;
    }

    /// The error follows n_max * max|value|^2, which the output bound n_min * max|v1| * max|v2| underestimates
    vector<long long> product(vector<long long>& v1, vector<long long>& v2, bool same){
        long long max_value = 0;
        for (auto x: v1) max_value = max(max_value, abs(x));
        for (auto x: v2) max_value = max(max_value, abs(x));

        long double bound = (long double)max(v1.size(), v2.size()) * max_value * max_value;
        if (bound <= DOUBLE_LIMIT) return real_product(fast, v1, v2, same);
        if (bound <= EXACT_LIMIT) return real_product(wide, v1, v2, same);
        assert(max_value < LL_MULTIPLY_LIMIT);
        if (split_fits_double(v1, v2)) return split_product(fast, v1, v2);
        return split_product(wide, v1, v2);
    }

    /// Real data packs into half-length complex z: E = (Z[k] + conj Z[-k]) / 2, O = (Z[k] - conj Z[-k]) / 2i, X[k] = E + w^k O, X[k + len / 2] = E - w^k O
    template <typename T>
    vector<long long> real_product(Engine<T>& e, vector<long long>& v1, vector<long long>& v2, bool same){
        int p_len = v1.size() + v2.size() - 1, len = pad(v1, v2), half = len / 2;
        e.prepare(half, len);

        for (int k = 0; k < half; k++) e.u[k] = Complex<T>(v1[2 * k], v1[2 * k + 1]);
        e.transform(e.u, e.f, half);
        if (!same){
            for (int k = 0; k < half; k++) e.w[k] = Complex<T>(v2[2 * k], v2[2 * k + 1]);
            e.transform(e.w, e.g, half);
        }
        vector<Complex<T>>& h = same ? e.f : e.g;

        const Complex<T> half_real(0.5, 0), half_imag(0, -0.5), unit_imag(0, 1);
        for (int k = 0; k < half; k++){
            int j = (half - k) & (half - 1);
            Complex<T> wk = e.roots[half + k];
            Complex<T> even1 = (e.f[k] + e.f[j].conjugate()) * half_real, odd1 = (e.f[k] - e.f[j].conjugate()) * half_imag * wk;
            Complex<T> even2 = (h[k] + h[j].conjugate()) * half_real, odd2 = (h[k] - h[j].conjugate()) * half_imag * wk;
            Complex<T> low = (even1 + odd1) * (even2 + odd2), high = (even1 - odd1) * (even2 - odd2);
            e.u[k] = (low + high) * half_real + unit_imag * ((low - high) * half_real * wk.conjugate());
        }
        e.transform(e.u, e.f, half);

        /// A forward transform read at index -m is the inverse transform scaled by half
        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            const Complex<T>& c = e.f[(half - (i >> 1)) & (half - 1)];
            res[i] = round_to_nearest(((i & 1) ? c.img : c.real) / half);
        }

        return res;
    }

    /// The split halves stay below 2^15, so their products need max(n, m) * 2^30 <= DOUBLE_LIMIT to run in double
    static bool split_fits_double(const vector<long long>& v1, const vector<long long>& v2){
        return (long double)max(v1.size(), v2.size()) * (1LL << (2 * SPLIT_BITS)) <= DOUBLE_LIMIT;
    }

    /// Negative values split into floor(v / 2^15) and v mod 2^15, and recombining modulo 2^64 is exact when the coefficient fits a long long
    template <typename T>
    vector<long long> split_product(Engine<T>& e, vector<long long>& v1, vector<long long>& v2){
        int p_len = v1.size() + v2.size() - 1, len = split_multiply(e, v1, v2);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            unsigned long long low = round_to_nearest(e.f[i].real), mid = round_to_nearest(e.g[i].real);
            unsigned long long high = round_to_nearest(e.f[i].img);
            res[i] = low + (mid << SPLIT_BITS) + (high << (2 * SPLIT_BITS));
        }

        return res;
    }

    template <typename T>
    vector<long long> split_mod(Engine<T>& e, vector<long long>& v1, vector<long long>& v2, int mod){
        int p_len = v1.size() + v2.size() - 1, len = split_multiply(e, v1, v2);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            long long low = round_to_nearest(e.f[i].real), mid = round_to_nearest(e.g[i].real), high = round_to_nearest(e.f[i].img);
            res[i] = (low + ((mid % mod) << SPLIT_BITS) + ((high % mod) << (2 * SPLIT_BITS))) % mod;
        }

        return res;
    }

    /// Leaves low * low, the cross terms and high * high in f.real, g.real and f.img, in 4 transforms or 3 for equal inputs
    template <typename T>
    int split_multiply(Engine<T>& e, vector<long long>& v1, vector<long long>& v2){
        const int mask = (1 << SPLIT_BITS) - 1;

        int len = pad(v1, v2);
        e.prepare(len, len);
        for (int i = 0; i < len; i++) e.u[i] = Complex<T>(v1[i] & mask, v1[i] >> SPLIT_BITS);
        for (int i = 0; i < len; i++) e.w[i] = Complex<T>(v2[i] & mask, v2[i] >> SPLIT_BITS);

        e.transform(e.u, e.f, len);
        if (v1 == v2){
            for (int i = 0; i < len; i++) e.g[i] = e.f[i];
        }
        else e.transform(e.w, e.g, len);

        for (int i = 0; i < len; i++){
            int j = (len - i) & (len - 1);
            Complex<T> c1 = e.f[j].conjugate(), c2 = e.g[j].conjugate();

            Complex<T> a1 = (e.f[i] + c1) * Complex<T>(0.5, 0);
            Complex<T> a2 = (e.f[i] - c1) * Complex<T>(0, -0.5);
            Complex<T> b1 = (e.g[i] + c2) * Complex<T>(0.5 / len, 0);
            Complex<T> b2 = (e.g[i] - c2) * Complex<T>(0, -0.5 / len);
            e.u[j] = a1 * b1 + a2 * b2 * Complex<T>(0, 1);
            e.w[j] = a1 * b2 + a2 * b1;
        }

        e.transform(e.u, e.f, len);
        e.transform(e.w, e.g, len);
        return len;
    }

    /// Pads v1 and v2 with zeros to the smallest power of two >= n + m, which is at least 2, so a half transform exists
    static int pad(vector<long long>& v1, vector<long long>& v2){
        int n = v1.size(), m = v2.size();
        while (n > 1 && v1[n - 1] == 0) n--;
        while (m > 1 && v2[m - 1] == 0) m--;

        int len = 1;
        while (len < n + m) len <<= 1;
        v1.resize(len, 0), v2.resize(len, 0);
        return len;
    }

    static long long round_to_nearest(long double x){
        long long res = fabsl(x) + 0.5L;
        return x < 0 ? -res : res;
    }

    static bool is_binary_string(const string& str){
        return all_of(str.begin(), str.end(), [](char c){ return c == '0' || c == '1'; });
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

#define MAX 400010
#define MOD 1000003

int n;
long long a, b, c, L[MAX], T[MAX], fact[MAX], inv[MAX], pa[MAX], pb[MAX];

long long expo(long long x, long long k){
    long long res = 1;
    for (x %= MOD; k; k >>= 1){
        if (k & 1) res = res * x % MOD;
        x = x * x % MOD;
    }
    return res;
}

long long binomial(int x, int k){
    return fact[x] * inv[k] % MOD * inv[x - k] % MOD;
}

int main(){
    using namespace fio;
    int i, k;

    read(n, a, b, c);
    for (i = 1; i <= n; i++) read(L[i]);
    for (i = 1; i <= n; i++) read(T[i]);

    fact[0] = 1;
    for (i = 1; i <= 2 * n; i++) fact[i] = fact[i - 1] * i % MOD;
    inv[2 * n] = expo(fact[2 * n], MOD - 2);
    for (i = 2 * n; i > 0; i--) inv[i - 1] = inv[i] * i % MOD;

    pa[0] = pb[0] = 1;
    for (i = 1; i <= n; i++) pa[i] = pa[i - 1] * a % MOD, pb[i] = pb[i - 1] * b % MOD;

    long long res = 0;
    for (k = 2; k <= n; k++){
        long long ways = binomial(2 * n - k - 2, n - 2);
        res += L[k] * pa[n - 1] % MOD * pb[n - k] % MOD * ways % MOD;
        res += T[k] * pb[n - 1] % MOD * pa[n - k] % MOD * ways % MOD;
    }

    vector<long long> u(n - 1), v(n - 1);
    for (i = 0; i <= n - 2; i++) u[i] = pa[i] * inv[i] % MOD, v[i] = pb[i] * inv[i] % MOD;

    FFT fft;
    vector<long long> w = fft.mod_multiply(u, v, MOD);
    for (i = 0; i < (int)w.size(); i++) res += c * fact[i] % MOD * w[i] % MOD;

    write(res % MOD);
    return 0;
}
