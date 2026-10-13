/***
 * Fast Fourier Transformation
 *
 * Uses custom class for complex numbers and various optimizations for better performance
 *
 * Usage (FFT fft; each method is documented above its definition):
 *   fft.multiply(a, b)               // exact product of small integer vectors, negatives allowed
 *   fft.square(a)                    // same as multiply(a, a) but faster
 *   fft.mod_multiply(a, b, mod)      // product modulo mod, 1 <= mod < 2^31, results in [0, mod)
 *   fft.ll_multiply(a, b)            // exact product of values in [0, LL_MULTIPLY_LIMIT)
 *   fft.hamming_distance(str, pat)   // mismatches of binary pat against each window of str, |str| - |pat| + 1 values
 *   fft.common_ones(str, pat)        // common 1 bits of binary pat and each window of str, |str| - |pat| + 1 values
 * All products are linear (no wrap around) of length |a| + |b| - 1, an empty input gives an empty result
 * The two string methods return an empty vector when |pat| > |str|
 *
 * Circular convolution of two size n vectors: fold the linear product, add % mod for mod_multiply
 *   auto c = fft.multiply(a, b);
 *   for (int i = n; i < (int)c.size(); i++) c[i - n] += c[i];
 *   c.resize(n);
 *
 * Complexity:
 *   - O(len) to extend the roots of unity the first time a transform length is used
 *   - O(n log n) for all exposed methods
 *
 * Memory: len is the smallest power of two >= n + m, an instance keeps its roots of unity and work buffers
 * at the largest len it has used, 2 complex arrays of len (4 once mod_multiply or ll_multiply ran)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Use double to gain more speed at the cost of precision (double should be fine for most problems)
typedef long double fType;

/// Work buffers are members reused across calls: allocating them per call measured 10-25% slower at len = 2^20
struct FFT{
    static constexpr int MOD_SPLIT_LIMIT = 15;
    static constexpr int LL_MULTIPLY_LIMIT = 1500000000;

    /***
     * Same as multiply(v, v) but faster
     *
    ***/
    vector<long long> square(vector<long long> v){
        if (v.empty()) return {};
        int i, len, p_len = 2 * v.size() - 1;

        len = build(v, v);
        for (i = 0; i < len; i++) u[i] = ComplexNum(v[i], 0);
        transform(u, f, len);
        for (i = 0; i < len; i++) u[i] = f[i] * f[i];
        transform(u, f, len);

        /// A forward transform read at index -i is the inverse transform scaled by len
        vector<long long> res(p_len, 0);
        for (i = 0; i < min(len, p_len); i++){
            res[i] = round_to_nearest(f[(len - i) & (len - 1)].real / (fType)len);
        }

        return res;
    }

    /***
     * Multiplies two polynomials and returns the result in a vector
     * Basically an optimized version of the following:
     *
     * vector multiply(vector v1, vector v2):
     *     n = v1.size()
     *     m = v2.size()
     *     vector res = vector(n + m - 1, 0)
     *
     *     for (i = 0; i < n; i++)
     *         for (j = 0; j < m; j++)
     *             res[i + j] += v1[i] * v2[j]
     *
     *     return res
     *
     * Note that the values in vector shouldn't be too large in general
     * Multiplying them and adding them up may cause precision error otherwise
     * The error margin will depend on two things mostly:
     *     i) The magnitude of the initial numbers
     *    ii) The chosen data type for complex numbers - double or long double, etc
     *
     * In general, if the values are larger than 10^6, it'd be better to use ll_multiply instead
     *
    ***/
    vector<long long> multiply(vector<long long> v1, vector<long long> v2){
        if (v1.empty() || v2.empty()) return {};
        if (is_equal(v1, v2)) return square(v1);

        int i, j, len, p_len = v1.size() + v2.size() - 1;

        len = build(v1, v2);
        for (i = 0; i < len; i++) u[i] = ComplexNum(v1[i], v2[i]);
        transform(u, f, len);

        for (i = 0; i < len; i++){
            j = (len - 1) & (len - i);
            u[i] = (f[j] * f[j] - f[i].conjugate() * f[i].conjugate()) * ComplexNum(0, -0.25 / len);
        }
        transform(u, f, len);

        vector<long long> res(p_len, 0);
        for (i = 0; i < min(len, p_len); i++){
            res[i] = round_to_nearest(f[i].real);
        }

        return res;
    }

    /***
     * Same as multiply(v1, v2), only values are calculated modulo mod, 1 <= mod < 2^31
     * Inputs may be negative or >= mod, they are reduced first, results lie in [0, mod)
     *
    ***/
    vector<long long> mod_multiply(vector<long long> v1, vector<long long> v2, int mod){
        if (v1.empty() || v2.empty()) return {};
        const int bits = MOD_SPLIT_LIMIT;
        int p_len = v1.size() + v2.size() - 1;
        for (auto&& x: v1) x = (x % mod + mod) % mod;
        for (auto&& x: v2) x = (x % mod + mod) % mod;

        int len = split_multiply(v1, v2, bits);
        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            long long x = round_to_nearest(f[i].real);
            long long y = round_to_nearest(g[i].real);
            long long z = round_to_nearest(f[i].img);
            res[i] = (x + ((y % mod) << bits) + ((z % mod) << (2 * bits))) % mod;
        }

        return res;
    }

    /***
     * Same as multiply(v1, v2), but supports larger values
     * Unlike multiply, does not support negative values as the values are split into 15 bit halves
     * Values in the vectors should not exceed the constant LL_MULTIPLY_LIMIT
     * Every coefficient of the result must be below 2^63
     *
    ***/
    vector<long long> ll_multiply(vector<long long> v1, vector<long long> v2){
        for (auto x: v1) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        for (auto x: v2) assert(x >= 0 && x < LL_MULTIPLY_LIMIT);
        if (v1.empty() || v2.empty()) return {};

        int p_len = v1.size() + v2.size() - 1;
        int len = split_multiply(v1, v2, MOD_SPLIT_LIMIT);

        vector<long long> res(p_len, 0);
        for (int i = 0; i < min(len, p_len); i++){
            unsigned long long x = round_to_nearest(f[i].real);
            unsigned long long y = round_to_nearest(g[i].real);
            unsigned long long z = round_to_nearest(f[i].img);
            res[i] = x + (y << MOD_SPLIT_LIMIT) + (z << (2 * MOD_SPLIT_LIMIT));
        }

        return res;
    }

    /***
     * Hamming distance vector with every substring of length |pattern| in str
     * Hamming distance is basically the number of mismatches
     *
     * str and pattern consists of only '1' and '0'
     *
     * For example,
     *     str: "1000100101"
     *     pattern: "0110"
     *     returns: [3, 3, 1, 1, 4, 1, 2]
     *
    ***/
    vector<long long> hamming_distance(const string& str, const string& pattern){
        assert(is_binary_string(str) && is_binary_string(pattern));

        int i, j, n = str.size(), m = pattern.size();
        vector<long long> res, v1(n, 0), v2(m, 0);
        if (m == 0) return vector<long long>(n + 1, 0);

        for (i = 0; i < n; i++) v1[i] = str[i] == '1' ? 1 : -1;
        for (i = 0, j = m - 1; j >= 0; i++, j--) v2[i] = pattern[j] == '1' ? 1 : -1;

        auto v = multiply(v1, v2);
        for (i = 0; (i + m) <= n; i++){
            res.push_back(m - ((v[i + m - 1] + m) >> 1));
        }

        return res;
    }

    /***
     * Number of common 1 bits of pattern and every substring of length |pattern| in str
     *
     * str and pattern consists of only '1' and '0'
     *
     * For example,
     *     str: "0110110"
     *     pattern: "110"
     *     returns: [1, 2, 1, 1, 2]
     *
    ***/
    vector<long long> common_ones(const string& str, const string& pattern){
        assert(is_binary_string(str) && is_binary_string(pattern));

        int i, n = str.size(), m = pattern.size();
        vector<long long> v1(n, 0), v2(m, 0);
        if (m == 0) return vector<long long>(n + 1, 0);

        for (i = 0; i < n; i++) v1[i] = str[n - i - 1] - 48;
        for (i = 0; i < m; i++) v2[i] = pattern[i] - 48;

        vector<long long> res;
        auto v = multiply(v1, v2);
        for (i = n; i >= m; i--) res.push_back(v[i - 1]);
        return res;
    }

    struct ComplexNum{
        fType real, img;

        ComplexNum() {}  /// left uninitialised so growing a work buffer does not write it twice, every entry is written before it is read
        ComplexNum(fType real, fType img=0) : real(real), img(img) {}

        inline ComplexNum conjugate(){
            return ComplexNum(real, -img);
        }

        inline ComplexNum operator + (ComplexNum x){
            return ComplexNum(real + x.real, img + x.img);
        }

        inline ComplexNum operator - (ComplexNum x){
            return ComplexNum(real - x.real, img - x.img);
        }

        inline ComplexNum operator * (ComplexNum x){
            return ComplexNum(real * x.real - img * x.img, real * x.img + img * x.real);
        }
    };

    int last_len = -1;
    vector<int> rev;
    vector<ComplexNum> u, w, f, g;

    /// roots[k + j] = e^(i * pi * j / k) for every power of two k < roots.size(), grown on demand
    vector<ComplexNum> roots = {ComplexNum(0), ComplexNum(1)};

    static long long round_to_nearest(const fType& x){
        long long res = abs(x) + 0.5;
        return (x < 0) ? -res : res;
    }

    static int get_bit(int len){
        return 32 - __builtin_clz(len) - (__builtin_popcount(len) == 1);
    }

    static bool is_equal(const vector<long long>& v1, const vector<long long>& v2){
        if (v1.size() != v2.size()) return false;

        for (int i = 0; i < (int)v1.size(); i++){
            if (v1[i] != v2[i]) return false;
        }

        return true;
    }

    static bool is_binary_string(const string& str){
        return all_of(str.begin(), str.end(), [](char c){ return c == '0' || c == '1'; });
    }

    void extend_roots(int len){
        roots.reserve(len);  /// one allocation instead of a copy per doubling
        for (int k = roots.size(); k < len; k <<= 1){
            fType theta = acosl(-1.0L) / k;
            ComplexNum mul = ComplexNum(cos(theta), sin(theta));

            roots.resize(2 * k);
            for (int j = k >> 1; j < k; j++){
                roots[2 * j] = roots[j];
                roots[2 * j + 1] = roots[j] * mul;
            }
        }
    }

    /// Pads v1 and v2 to the transform length len, grows u, f and rev to len and returns len
    int build(vector<long long>& v1, vector<long long>& v2){
        int i, n = v1.size(), m = v2.size();
        while (n > 1 && v1[n - 1] == 0) n--;
        while (m > 1 && v2[m - 1] == 0) m--;

        int len = 1 << get_bit(n + m);
        v1.resize(len, 0), v2.resize(len, 0);

        extend_roots(len);
        if ((int)u.size() < len) u.resize(len), f.resize(len), rev.resize(len);
        if (len != last_len){
            last_len = len;
            const int bit = get_bit(len);
            for (i = 1; i < len; i++){
                rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (bit - 1));
            }
        }

        return len;
    }

    void transform(const vector<ComplexNum>& in, vector<ComplexNum>& out, int len){
        for (int i = 0; i < len; i++) out[i] = in[rev[i]];

        for (int k = 1; k < len; k <<= 1){
            for (int i = 0; i < len; i += (k << 1)){
                ComplexNum z, *a = &out[i], *b = &out[i + k], *c = &roots[k];
                if (k == 1){
                    z = (*b) * (*c);
                    *b = *a - z, *a = *a + z;
                }

                for (int j = 0; j < k && k > 1; j += 2, a++, b++, c++){
                    z = (*b) * (*c);
                    *b = *a - z, *a = *a + z;
                    a++, b++, c++;
                    z = (*b) * (*c);
                    *b = *a - z, *a = *a + z;
                }
            }
        }
    }

    /// Splits values into low and high bits, leaving low * low, the cross terms and high * high
    /// in f.real, g.real and f.img respectively, 4 transforms in total (3 for equal inputs)
    int split_multiply(vector<long long>& v1, vector<long long>& v2, int bits){
        const int mask = (1 << bits) - 1;

        int i, j, len = build(v1, v2);
        if ((int)w.size() < len) w.resize(len), g.resize(len);
        for (i = 0; i < len; i++) u[i] = ComplexNum(v1[i] & mask, v1[i] >> bits);
        for (i = 0; i < len; i++) w[i] = ComplexNum(v2[i] & mask, v2[i] >> bits);

        transform(u, f, len);
        for (i = 0; i < len; i++) g[i] = f[i];
        if (!is_equal(v1, v2)) transform(w, g, len);

        for (i = 0; i < len; i++){
            j = (len - 1) & (len - i);
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
};

int main(){
    FFT fft;

    vector<long long> v1, v2, expected_result;

    v1 = {5, 1, 2, 6, 9, 8};
    v2 = {3, 9, 0, 2};
    assert(fft.multiply(v1, v2) == vector<long long>({15, 48, 15, 46, 83, 109, 84, 18, 16}));

    int mod = 14;
    assert(fft.mod_multiply(v1, v2, mod) == vector<long long>({1, 6, 1, 4, 13, 11, 0, 4, 2}));

    for (auto && x: v1) x = (1 << 30) - x;
    for (auto && x: v2) x = (1 << 30) - x;
    expected_result = {1152921496016912399, 2305842989886341168, 3458764492345704463, 4611685988362616878, 4611685984067649619, 4611685976551456877, 3458764477313318996, 2305842988812599314, 1152921493869428752};

    assert(fft.ll_multiply(v1, v2) == expected_result);

    assert(fft.mod_multiply({-1, 2}, {1, 1}, 7) == vector<long long>({6, 1, 2}));
    assert(fft.multiply({}, {1, 2}).empty() && fft.mod_multiply({}, {}, 7).empty());

    expected_result = {3, 3, 1, 1, 4, 1, 2};
    assert(fft.hamming_distance("1000100101", "0110") == expected_result);

    expected_result = {1, 2, 1, 1, 2};
    assert(fft.common_ones(string("0110110"), "110") == expected_result);

    expected_result = {0, 0, 0, 0, 0};
    assert(fft.hamming_distance("0101", "") == expected_result && fft.common_ones("0101", "") == expected_result);
    assert(fft.hamming_distance("", "") == vector<long long>({0}) && fft.common_ones("01", "011").empty());

    return 0;
}
