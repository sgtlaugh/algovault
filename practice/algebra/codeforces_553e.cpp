/***
 *
 * Fast Fourier Transform
 * Polynomial multiplication, exact over integers or modulo any mod, and binary string matching
 *
 * Complexity:
 *   - every call: O((n + m) log(n + m)) for inputs of sizes n and m
 *   - memory: O(largest n + m so far), kept between calls
 *
 * FFT fft: one instance for any number of calls
 * fft.common_ones(str, pat): common 1 bits of binary pat and each window of str, |str| - |pat| + 1 values
 * fft.hamming_distance(str, pat): mismatches of binary pat against each window of str, |str| - |pat| + 1 values
 * fft.ll_multiply(a, b): the exact product of values in [0, LL_MULTIPLY_LIMIT), every coefficient below 2^63
 * fft.mod_multiply(a, b, mod): the product modulo mod, 1 <= mod < 2^31, any input values, results in [0, mod)
 * fft.multiply(a, b): the product, |values| < LL_MULTIPLY_LIMIT, exact while every coefficient fits in a long long
 * fft.square(a): the same as multiply(a, a)
 *
 * Notes:
 *   - products have length n + m - 1, and an empty input gives an empty result
 *   - the string methods return an empty vector when |pat| > |str| and |str| + 1 zeros for an empty pat
 *   - use ntt.cpp modulo an NTT-friendly prime like 998244353 for improved speed
 *   - mod_multiply with mod near 10^9 and an input longer than about 1.8e5 uses long double, ntt.cpp is about 2x faster
 *   - transforms round n + m - 1 up to a power of two, so one coefficient past it doubles the work: drop unneeded terms
 *   - for circular convolution of two size n vectors: fold the linear product, add % mod for mod_multiply
 *         auto c = fft.multiply(a, b);
 *         for (int i = n; i < (int)c.size(); i++) c[i - n] += c[i];
 *         c.resize(n);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FFT{
    template <typename T>
    struct Complex{
        T real, img;

        Complex() {}
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

        vector<Complex<T>> u, w, f, g;
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
            /// The first stage's twiddle is 1, so it is done while gathering in bit-reversed order
            if (len == 1) out[0] = in[0];
            for (int i = 0; i + 1 < len; i += 2){
                Complex<T> x = in[rev[i]], y = in[rev[i + 1]];
                out[i] = x + y, out[i + 1] = x - y;
            }

            for (int k = 2; k < len; k <<= 1){
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

    /// Largest n_max * max|value|^2 one transform keeps exact, at least 4x below the first error up to n = 10^6
    static constexpr double DOUBLE_LIMIT = 0.01 / numeric_limits<double>::epsilon();
    static constexpr long double EXACT_LIMIT = 0.01L / numeric_limits<long double>::epsilon();

    Engine<double> fast;
    Engine<long double> wide;

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

    vector<long long> ll_multiply(vector<long long> v1, vector<long long> v2){
        for (auto x: v1) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        for (auto x: v2) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        if (v1.empty() || v2.empty()) return {};
        return product(v1, v2, v1 == v2);
    }

    vector<long long> mod_multiply(vector<long long> v1, vector<long long> v2, int mod){
        if (v1.empty() || v2.empty()) return {};
        /// Residues in (-mod / 2, mod / 2] and balanced digits halve every magnitude, a 4x smaller rounding error
        for (auto& x: v1) x = centered(x, mod);
        for (auto& x: v2) x = centered(x, mod);

        if ((long double)max(v1.size(), v2.size()) * (mod / 2) * (mod / 2) <= EXACT_LIMIT){
            vector<long long> res = product(v1, v2, v1 == v2);
            for (auto& x: res) x = non_negative(x % mod, mod);
            return res;
        }

        long long base = sqrtl(mod);
        while (base * base < mod) base++;
        if (split_fits_double(v1, v2, base / 2 + 1)) return split_mod(fast, v1, v2, mod, base);
        return split_mod(wide, v1, v2, mod, base);
    }

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

    static bool is_binary_string(const string& str){
        return all_of(str.begin(), str.end(), [](char c){ return c == '0' || c == '1'; });
    }

    vector<long long> product(vector<long long>& v1, vector<long long>& v2, bool same){
        long long max_value = 0;

        for (auto x: v1) max_value = max(max_value, abs(x));
        for (auto x: v2) max_value = max(max_value, abs(x));

        long double bound = (long double)max(v1.size(), v2.size()) * max_value * max_value;
        if (bound <= DOUBLE_LIMIT) return real_product(fast, v1, v2, same);
        if (bound <= EXACT_LIMIT) return real_product(wide, v1, v2, same);
        assert(max_value < LL_MULTIPLY_LIMIT);

        long long max_part = max(1LL << SPLIT_BITS, (max_value >> SPLIT_BITS) + 1);
        if (split_fits_double(v1, v2, max_part)) return split_product(fast, v1, v2);
        return split_product(wide, v1, v2);
    }

    static bool split_fits_double(const vector<long long>& v1, const vector<long long>& v2, long long max_part){
        return (long double)max(v1.size(), v2.size()) * max_part * max_part <= DOUBLE_LIMIT;
    }

    template <typename T>
    vector<long long> split_mod(Engine<T>& e, vector<long long>& v1, vector<long long>& v2, int mod, long long base){
        int p_len = v1.size() + v2.size() - 1, len = split_multiply(e, v1, v2, base);
        long long base_sq = base * base % mod;

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            long long low = round_to_nearest(e.f[i].real), mid = round_to_nearest(e.g[i].real), high = round_to_nearest(e.f[i].img);
            /// below 4.7e18 as high % mod * base_sq < 2^62
            res[i] = non_negative((low + mid % mod * base + high % mod * base_sq) % mod, mod);
        }

        return res;
    }

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

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            const Complex<T>& c = e.f[(half - (i >> 1)) & (half - 1)];
            res[i] = round_to_nearest(((i & 1) ? c.img : c.real) / half);
        }

        return res;
    }

    template <typename T>
    vector<long long> split_product(Engine<T>& e, vector<long long>& v1, vector<long long>& v2){
        int p_len = v1.size() + v2.size() - 1, len = split_multiply(e, v1, v2, 1LL << SPLIT_BITS);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            unsigned long long low = round_to_nearest(e.f[i].real), mid = round_to_nearest(e.g[i].real);
            unsigned long long high = round_to_nearest(e.f[i].img);
            res[i] = low + (mid << SPLIT_BITS) + (high << (2 * SPLIT_BITS));
        }

        return res;
    }

    template <typename T>
    int split_multiply(Engine<T>& e, vector<long long>& v1, vector<long long>& v2, long long base){
        int len = pad(v1, v2);

        e.prepare(len, len);
        for (int i = 0; i < len; i++) e.u[i] = split<T>(v1[i], base);
        for (int i = 0; i < len; i++) e.w[i] = split<T>(v2[i], base);

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

    static long long centered(long long x, int mod){
        if (x < 0 || x >= mod) x = (x % mod + mod) % mod;
        return x > mod / 2 ? x - mod : x;
    }

    static long long non_negative(long long x, int mod){
        return x < 0 ? x + mod : x;
    }

    template <typename T>
    static long long round_to_nearest(T x){
        long long res = abs(x) + T(0.5);
        return x < 0 ? -res : res;
    }

    static int pad(vector<long long>& v1, vector<long long>& v2){
        int n = v1.size(), m = v2.size();

        while (n > 1 && v1[n - 1] == 0) n--;
        while (m > 1 && v2[m - 1] == 0) m--;

        int len = 1;
        while (len < n + m) len <<= 1;
        v1.resize(len, 0), v2.resize(len, 0);
        return len;
    }

    /// Split values are below LL_MULTIPLY_LIMIT or mod < 2^31 in absolute value, so cheaper 32-bit division is exact
    template <typename T>
    static Complex<T> split(int32_t v, int32_t base){
        int32_t low = v % base, high = v / base, carry = (2 * low > base) - (2 * low < -base);
        return Complex<T>(low - carry * base, high + carry);
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

#define MAXT 20010
#define MAXM 110
#define SHIFT 30
#define LIMB 28
#define BLOCK 128

FFT fft;
int n, m, t, x, from[MAXM], to[MAXM], cost[MAXM], p[MAXM][MAXT];
long long dist[60], grid[60][MAXT];
long double tail[MAXM][MAXT], best[60][MAXT];
__int128 acc[MAXM][MAXT];

void leaf(int tau){
    int v, e;

    for (v = 1; v <= n; v++) best[v][tau] = (v == n) ? 0 : 1e18;
    for (e = 0; e < m; e++){
        if (from[e] == n) continue;
        long double cur = cost[e] + tail[e][tau] + (long double)acc[e][tau] / (100000.0L * (1 << SHIFT));
        best[from[e]][tau] = min(best[from[e]][tau], cur);
    }
    for (v = 1; v <= n; v++) grid[v][tau] = llroundl(best[v][tau] * (1 << SHIFT));
}

void solve(int l, int r){
    int e, i, u, mid, len;

    if ((r - l) < BLOCK){
        for (i = r; i >= l; i--){
            leaf(i);
            for (e = 0; e < m; e++){
                for (u = l; u < i; u++) acc[e][u] += (__int128)p[e][i - u] * grid[to[e]][i];
            }
        }

        return;
    }

    mid = (l + r) >> 1, len = r - mid;
    solve(mid + 1, r);

    for (e = 0; e < m; e++){
        vector<long long> hi(len), lo(len), q(p[e], p[e] + (r - l + 1));
        for (i = 0; i < len; i++){
            long long g = grid[to[e]][r - i];
            hi[i] = g >> LIMB, lo[i] = g & ((1LL << LIMB) - 1);
        }

        vector<long long> c_hi = fft.multiply(hi, q), c_lo = fft.multiply(lo, q);
        for (u = 0; u <= mid - l; u++) acc[e][mid - u] += ((__int128)c_hi[len + u] << LIMB) + c_lo[len + u];
    }

    solve(l, mid);
}

int main(){
    int i, j, k, e;
    long long suffix;

    read(n, m, t, x);
    for (i = 1; i <= n; i++) dist[i] = (i == n) ? 0 : LLONG_MAX / 4;
    for (e = 0; e < m; e++){
        read(from[e], to[e], cost[e]);
        for (k = 1; k <= t; k++) read(p[e][k]);
    }

    for (i = 0; i < n; i++){
        for (e = 0; e < m; e++) dist[from[e]] = min(dist[from[e]], dist[to[e]] + cost[e]);
    }

    for (e = 0; e < m; e++){
        suffix = 0;
        for (j = 0; j <= t; j++){
            if (j) suffix += p[e][t - j + 1];
            tail[e][j] = (long double)suffix * (dist[to[e]] + x) / 100000.0L;
        }
    }

    solve(0, t);
    write(fixed(best[1][0], 10));
    return 0;
}
