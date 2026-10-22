/***
 *
 * Fast Fourier Transform
 * Polynomial multiplication, exact over integers or modulo any mod, and binary string matching
 *
 * Complexity:
 *   - every call: O((n + m) log(n + m)) for inputs of sizes n and m
 *   - the first call at a new transform length len: O(len) more to extend the roots of unity
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
 * multiply and square use one transform while max(n, m) * max|value|^2 <= EXACT_LIMIT, about 9.2e16 with long double
 * Past EXACT_LIMIT they split values into 15 bit halves like ll_multiply, about 2x slower, and need |values| < LL_MULTIPLY_LIMIT
 * The string methods return an empty vector when |pat| > |str| and |str| + 1 zeros for an empty pat
 * Memory: 2 complex arrays of the largest len used, the smallest power of two >= n + m, and 4 once a split product ran
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
    /// double is faster, but one transform is then exact only to about 1e14 instead of 1e17 (see EXACT_LIMIT)
    typedef long double fType;

    struct ComplexNum{
        fType real, img;

        ComplexNum() {}  /// left uninitialised so growing a work buffer does not write it twice, every entry is written before it is read
        ComplexNum(fType real, fType img=0) : real(real), img(img) {}

        ComplexNum conjugate(){
            return ComplexNum(real, -img);
        }

        ComplexNum operator + (ComplexNum x){
            return ComplexNum(real + x.real, img + x.img);
        }

        ComplexNum operator - (ComplexNum x){
            return ComplexNum(real - x.real, img - x.img);
        }

        ComplexNum operator * (ComplexNum x){
            return ComplexNum(real * x.real - img * x.img, real * x.img + img * x.real);
        }
    };

    static constexpr int SPLIT_BITS = 15;
    static constexpr int LL_MULTIPLY_LIMIT = 1500000000;

    /// Largest n_max * max|value|^2 multiply and square still do with one transform, about 9.2e16 for long double and 4.5e13 for double
    /// Measured first errors up to n = 10^6: 6.3e17 for long double and 4e14 for double, so the limit keeps at least a 6x margin
    static constexpr fType EXACT_LIMIT = 0.01 / numeric_limits<fType>::epsilon();

    int last_len = -1;
    vector<int> rev;

    /// Work buffers are members reused across calls: allocating them per call measured 10-25% slower at len = 2^20
    vector<ComplexNum> u, w, f, g;

    /// roots[k + j] = e^(i * pi * j / k) for every power of two k < roots.size(), grown on demand
    vector<ComplexNum> roots = {ComplexNum(0), ComplexNum(1)};

    vector<long long> multiply(vector<long long> v1, vector<long long> v2){
        if (v1.empty() || v2.empty()) return {};
        if (v1 == v2) return square(v1);
        if (!fits_fast_path(v1, v2)) return split_product(v1, v2);

        int p_len = v1.size() + v2.size() - 1, len = build(v1, v2);
        for (int i = 0; i < len; i++) u[i] = ComplexNum(v1[i], v2[i]);
        transform(u, f, len);

        /// With f = FFT(v1 + i v2), FFT(v1 * v2)[j] = (f[j]^2 - conj(f[-j])^2) / 4i, stored at -j so a forward transform inverts it
        for (int i = 0; i < len; i++){
            int j = (len - i) & (len - 1);
            u[i] = (f[j] * f[j] - f[i].conjugate() * f[i].conjugate()) * ComplexNum(0, -0.25 / len);
        }
        transform(u, f, len);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++) res[i] = round_to_nearest(f[i].real);
        return res;
    }

    vector<long long> square(vector<long long> v){
        if (v.empty()) return {};
        if (!fits_fast_path(v, v)) return split_product(v, v);

        int p_len = 2 * v.size() - 1, len = build(v, v);
        for (int i = 0; i < len; i++) u[i] = ComplexNum(v[i], 0);
        transform(u, f, len);
        for (int i = 0; i < len; i++) u[i] = f[i] * f[i];
        transform(u, f, len);

        /// A forward transform read at index -i is the inverse transform scaled by len
        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++) res[i] = round_to_nearest(f[(len - i) & (len - 1)].real / (fType)len);
        return res;
    }

    vector<long long> mod_multiply(vector<long long> v1, vector<long long> v2, int mod){
        if (v1.empty() || v2.empty()) return {};
        for (auto& x: v1) x = (x % mod + mod) % mod;
        for (auto& x: v2) x = (x % mod + mod) % mod;

        int p_len = v1.size() + v2.size() - 1, len = split_multiply(v1, v2);
        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            long long low = round_to_nearest(f[i].real), mid = round_to_nearest(g[i].real), high = round_to_nearest(f[i].img);
            res[i] = (low + ((mid % mod) << SPLIT_BITS) + ((high % mod) << (2 * SPLIT_BITS))) % mod;
        }

        return res;
    }

    vector<long long> ll_multiply(vector<long long> v1, vector<long long> v2){
        for (auto x: v1) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        for (auto x: v2) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        if (v1.empty() || v2.empty()) return {};
        return split_product(v1, v2);
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

    /// multiply packs v1 + i v2 into one transform, so its error follows the larger input: n_max * max(|v1|, |v2|)^2
    /// The output bound n_min * max|v1| * max|v2| is not enough: n = 10^5, v1 in {0, 1}, v2 < 1.3 * 10^7 already errs
    static bool fits_fast_path(const vector<long long>& v1, const vector<long long>& v2){
        long long max_value = 0;
        for (auto x: v1) max_value = max(max_value, abs(x));
        for (auto x: v2) max_value = max(max_value, abs(x));
        if ((fType)max(v1.size(), v2.size()) * max_value * max_value <= EXACT_LIMIT) return true;

        assert(max_value < LL_MULTIPLY_LIMIT);
        return false;
    }

    /// A negative value splits into floor(v / 2^15) and v mod 2^15, so the parts stay as small as for |v|
    /// The parts are recombined modulo 2^64, exact whenever the true coefficient fits in a long long
    vector<long long> split_product(vector<long long>& v1, vector<long long>& v2){
        int p_len = v1.size() + v2.size() - 1, len = split_multiply(v1, v2);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            unsigned long long low = round_to_nearest(f[i].real), mid = round_to_nearest(g[i].real), high = round_to_nearest(f[i].img);
            res[i] = low + (mid << SPLIT_BITS) + (high << (2 * SPLIT_BITS));
        }

        return res;
    }

    /// Splits values into low and high bits, leaving low * low, the cross terms and high * high
    /// in f.real, g.real and f.img respectively, 4 transforms in total (3 for equal inputs)
    int split_multiply(vector<long long>& v1, vector<long long>& v2){
        const int mask = (1 << SPLIT_BITS) - 1;

        int len = build(v1, v2);
        if ((int)w.size() < len) w.resize(len), g.resize(len);
        for (int i = 0; i < len; i++) u[i] = ComplexNum(v1[i] & mask, v1[i] >> SPLIT_BITS);
        for (int i = 0; i < len; i++) w[i] = ComplexNum(v2[i] & mask, v2[i] >> SPLIT_BITS);

        transform(u, f, len);
        if (v1 == v2){
            for (int i = 0; i < len; i++) g[i] = f[i];
        }
        else transform(w, g, len);

        for (int i = 0; i < len; i++){
            int j = (len - i) & (len - 1);
            ComplexNum c1 = f[j].conjugate(), c2 = g[j].conjugate();

            ComplexNum a1 = (f[i] + c1) * ComplexNum(0.5, 0);
            ComplexNum a2 = (f[i] - c1) * ComplexNum(0, -0.5);
            ComplexNum b1 = (g[i] + c2) * ComplexNum(0.5 / len, 0);
            ComplexNum b2 = (g[i] - c2) * ComplexNum(0, -0.5 / len);
            u[j] = a1 * b1 + a2 * b2 * ComplexNum(0, 1);
            w[j] = a1 * b2 + a2 * b1;
        }

        transform(u, f, len);
        transform(w, g, len);
        return len;
    }

    /// Pads v1 and v2 to the transform length len, grows u, f and rev to len and returns len
    int build(vector<long long>& v1, vector<long long>& v2){
        int n = v1.size(), m = v2.size();
        while (n > 1 && v1[n - 1] == 0) n--;
        while (m > 1 && v2[m - 1] == 0) m--;

        int len = 1;
        while (len < n + m) len <<= 1;
        v1.resize(len, 0), v2.resize(len, 0);

        extend_roots(len);
        if ((int)u.size() < len) u.resize(len), f.resize(len), rev.resize(len);
        if (len != last_len){
            last_len = len;
            int bits = __lg(len);
            for (int i = 1; i < len; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bits - 1));
        }

        return len;
    }

    void extend_roots(int len){
        roots.reserve(len);
        for (int k = roots.size(); k < len; k <<= 1){
            fType theta = acosl(-1.0L) / k;
            ComplexNum mul(cos(theta), sin(theta));

            roots.resize(2 * k);
            for (int j = k >> 1; j < k; j++){
                roots[2 * j] = roots[j];
                roots[2 * j + 1] = roots[j] * mul;
            }
        }
    }

    void transform(const vector<ComplexNum>& in, vector<ComplexNum>& out, int len){
        for (int i = 0; i < len; i++) out[i] = in[rev[i]];

        for (int k = 1; k < len; k <<= 1){
            for (int i = 0; i < len; i += 2 * k){
                for (int j = 0; j < k; j++){
                    ComplexNum z = out[i + j + k] * roots[k + j];
                    out[i + j + k] = out[i + j] - z;
                    out[i + j] = out[i + j] + z;
                }
            }
        }
    }

    static long long round_to_nearest(fType x){
        long long res = abs(x) + 0.5;
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

#define MAXV 50000

int main(){
    using namespace fio;
    int a, b, c, i, j, s, v;
    string card, suits = "SHCD";
    vector<bool> composite(MAXV + 1);
    FFT fft;

    for (i = 2; (i * i) <= MAXV; i++){
        for (j = i * i; j <= MAXV; j += i) composite[j] = true;
    }

    while (read(a, b, c) == 3 && (a || b || c)){
        vector<vector<long long>> cards(4, vector<long long>(b + 1));
        for (s = 0; s < 4; s++){
            for (v = 0; v <= b; v++) cards[s][v] = composite[v];
        }

        while (c--){
            read(card);
            v = atoi(card.c_str());
            if (v <= b) cards[suits.find(card.back())][v] = 0;
        }

        vector<long long> ways = cards[0];
        for (s = 1; s < 4; s++){
            ways = fft.ll_multiply(ways, cards[s]);
            ways.resize(b + 1);
        }

        for (i = a; i <= b; i++) write(ways[i]);
        write("");
    }
    return 0;
}
