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

int main(){
    /// (5 + x + 2x^2 + 6x^3 + 9x^4 + 8x^5) * (3 + 9x + 2x^3)
    FFT fft;
    vector<long long> a = {5, 1, 2, 6, 9, 8}, b = {3, 9, 0, 2};
    assert(fft.multiply(a, b) == vector<long long>({15, 48, 15, 46, 83, 109, 84, 18, 16}));
    assert(fft.square({1, 2, 3}) == vector<long long>({1, 4, 10, 12, 9}));

    /// 2 * (10^9)^2 is past EXACT_LIMIT, so this takes the exact split product
    assert(fft.multiply({-1000000000, 3}, {1000000000, 1}) == vector<long long>({-1000000000000000000, 2000000000, 3}));

    /// Inputs are reduced first, -1 becomes 6 modulo 7
    assert(fft.mod_multiply(a, b, 14) == vector<long long>({1, 6, 1, 4, 13, 11, 0, 4, 2}));
    assert(fft.mod_multiply({-1, 2}, {1, 1}, 7) == vector<long long>({6, 1, 2}));

    assert(fft.ll_multiply({1000000000, 1}, {1000000000, 2}) == vector<long long>({1000000000000000000, 3000000000, 2}));

    /// Windows of "1000100101": 1000, 0001, 0010, 0100, 1001, 0010, 0101
    assert(fft.hamming_distance("1000100101", "0110") == vector<long long>({3, 3, 1, 1, 4, 1, 2}));
    assert(fft.common_ones("0110110", "110") == vector<long long>({1, 2, 1, 1, 2}));
    assert(fft.common_ones("01", "011").empty());

    return 0;
}
