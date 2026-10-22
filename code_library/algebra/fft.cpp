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
