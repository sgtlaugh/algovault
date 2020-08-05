/***
 * Fast Fourier Transformation
 *
 * Uses custom class for complex numbers and various optimizations for better performance
 *
 * Complexity: O(n log n) for all methods computing multiplication/convolution
 *
 * Notes:
 * Replace long double with double to gain more speed at the cost of precision
 * Double should be fine for most problems
 *
***/

#include <bits/stdc++.h>

#define MAX 2097152

using namespace std;

namespace fft{
    namespace{ /// Anonymous namespace to wrap internally used variables and methods
        const int MOD_SPLIT_LIMIT = 15;
        const int LL_MULTIPLY_MOD = 1500000000;

        long long A[MAX], B[MAX];
        int len, last = -1, step = 0, rev[MAX];

        struct ComplexNum{
            long double real, img;

            inline ComplexNum(long double real=0, long double img=0) : real(real), img(img) {}

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
        } u[MAX], v[MAX], f[MAX], g[MAX], dp[MAX], inv[MAX];

        inline long long round_to_nearest(long double x){
            long long res = abs(x) + 0.5;
            return (x < 0) ? -res : res;
        }

        void build(int& a, long long* A, int& b, long long* B){
            while (a > 1 && A[a - 1] == 0) a--;
            while (b > 1 && B[b - 1] == 0) b--;

            len = 1 << (32 - __builtin_clz(a + b) - (__builtin_popcount(a + b) == 1));
            for (int i = a; i < len; i++) A[i] = 0;
            for (int i = b; i < len; i++) B[i] = 0;

            if (!step++){
                dp[1] = inv[1] = ComplexNum(1);
                for (int i = 1; (1 << i) < MAX; i++){
                    long double theta = 2.0 * acos(0.0) / (1 << i);
                    ComplexNum mul = ComplexNum(cos(theta), sin(theta));
                    ComplexNum inv_mul = ComplexNum(cos(-theta), sin(-theta));

                    int lim = 1 << i;
                    for (int j = lim >> 1; j < lim; j++){
                        dp[2 * j] = dp[j], inv[2 * j] = inv[j];
                        inv[2 * j + 1] = inv[j] * inv_mul;
                        dp[2 * j + 1] = dp[j] * mul;
                    }
                }
            }

            if (last != len){
                last = len;
                int bit = (32 - __builtin_clz(len) - (__builtin_popcount(len) == 1));
                for (int i = 0; i < len; i++) rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (bit - 1));
            }
        }

        void transform(ComplexNum *in, ComplexNum *out, ComplexNum* ar){
            for (int i = 0; i < len; i++) out[i] = in[rev[i]];
            for (int k = 1; k < len; k <<= 1){
                for (int i = 0; i < len; i += (k << 1)){
                    ComplexNum z, *a = out + i, *b = out + i + k, *c = ar + k;
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

        int copy_to_array(long long *ar, const vector <long long>& v){
            int i, len = v.size();
            for (i = 0; i < len; i++) ar[i] = v[i];
            return len;
        }

        bool is_equal(const vector <long long>& v1, const vector <long long>& v2){
            if (v1.size() != v2.size()) return false;

            for (int i = 0; i < v1.size(); i++){
                if (v1[i] != v2[i]) return false;
            }
            return true;
        }

        bool is_binary_string(const char* str){
            for (int j = 0; str[j] != 0; j++){
                if (!(str[j] == '0' || str[j] == '1')) return false;
            }
            return true;
        }

        pair<vector<long long>, vector<long long>> build_convolution(const vector <long long>& v1, const vector <long long>& v2){
            assert(v1.size() == v2.size());

            int i, m, d = 0, n = v1.size();
            for (i = 0; i < n; i++) B[i] = B[i + n] = v2[i];
            for (i = 0; i < n; i++) A[i] = v1[i], A[i + n] = 0;
            n = 2 * n, m = 1 << (32 - __builtin_clz(n) - (__builtin_popcount(n) == 1));
            for (i = n; i < m; i++) A[i] = B[i] = 0;

            return {vector<long long>(A, A + n), vector<long long>(B, B + n)};
        }
    }

    /***
     * Same as multiply(v, v) but faster
     *
    ***/
    vector <long long> square(const vector <long long>& v){
        int a = copy_to_array(A, v);
        int p_len = 2 * a - 1;

        build(a, A, a, A);
        for (int i = 0; i < len; i++) u[i] = ComplexNum(A[i], 0);
        transform(u, f, dp);
        for (int i = 0; i < len; i++) u[i] = f[i] * f[i];
        transform(u, f, inv);

        for (int i = 0; i < len; i++) A[i] = round_to_nearest(f[i].real / (long double)len);
        for (int i = len; i < p_len; i++) A[i] = 0;

        return vector <long long>(A, A + p_len);
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
    vector <long long> multiply(const vector <long long>& v1, const vector <long long>& v2){
        if (is_equal(v1, v2)) return square(v1);

        int a = copy_to_array(A, v1);
        int b = copy_to_array(B, v2);
        int p_len = a + b - 1;

        build(a, A, b, B);
        for (int i = 0; i < len; i++) u[i] = ComplexNum(A[i], B[i]);
        transform(u, f, dp);

        for (int i = 0; i < len; i++){
            int j = (len - 1) & (len - i);
            u[i] = (f[j] * f[j] - f[i].conjugate() * f[i].conjugate()) * ComplexNum(0, -0.25 / len);
        }
        transform(u, f, dp);

        for (int i = 0; i < len; i++) A[i] = round_to_nearest(f[i].real);
        for (int i = len; i < p_len; i++) A[i] = 0;

        return vector <long long>(A, A + p_len);
    }

    /***
     * Same as multiply(v1, v2), only values are calculated modulo mod
     *
    ***/
    vector <long long> mod_multiply(const vector <long long>& v1, const vector <long long>& v2, int mod, int bits=MOD_SPLIT_LIMIT){
        const int mask = (1 << bits) - 1;

        int a = copy_to_array(A, v1);
        int b = copy_to_array(B, v2);
        int p_len = a + b - 1;

        build(a, A, b, B);
        for (int i = 0; i < len; i++) A[i] %= mod, B[i] %= mod;
        for (int i = 0; i < len; i++) u[i] = ComplexNum(A[i] & mask, A[i] >> bits);
        for (int i = 0; i < len; i++) v[i] = ComplexNum(B[i] & mask, B[i] >> bits);

        transform(u, f, dp);
        for (int i = 0; i < len; i++) g[i] = f[i];
        if (!is_equal(v1, v2)) transform(v, g, dp);

        for (int i = 0; i < len; i++){
            int j = (len - 1) & (len - i);
            ComplexNum c1 = f[j].conjugate(), c2 = g[j].conjugate();

            ComplexNum a1 = (f[i] + c1) * ComplexNum(0.5, 0);
            ComplexNum a2 = (f[i] - c1) * ComplexNum(0, -0.5);
            ComplexNum b1 = (g[i] + c2) * ComplexNum(0.5 / len, 0);
            ComplexNum b2 = (g[i] - c2) * ComplexNum(0, -0.5 / len);
            v[j] = a1 * b2 + a2 * b1;
            u[j] = a1 * b1 + a2 * b2 * ComplexNum(0, 1);
        }
        transform(u, f, dp);
        transform(v, g, dp);

        for (int i = 0; i < len; i++){
            long long x = round_to_nearest(f[i].real);
            long long y = round_to_nearest(g[i].real);
            long long z = round_to_nearest(f[i].img);
            A[i] = (x + ((y % mod) << bits) + ((z % mod) << (2 * bits))) % mod;
        }
        for (int i = len; i < p_len; i++) A[i] = 0;

        return vector <long long>(A, A + p_len);
    }

    /***
     * Same as multiply(v1, v2), but supports larger values
     * Unlike multiply, does not support negative values as it internally uses mod_multiply
     * Values in the vectors should not exceed the constant LL_MULTIPLY_MOD
     *
    ***/
    vector<long long> ll_multiply(const vector <long long>& v1, const vector <long long>& v2){
        for (auto x: v1) assert(x >= 0 && x < LL_MULTIPLY_MOD);
        for (auto x: v2) assert(x >= 0 && x < LL_MULTIPLY_MOD);

        int len = v1.size() + v2.size() - 1;
        const long long mod1 = LL_MULTIPLY_MOD, mod2 = mod1 + 1;

        auto mod_1_res = mod_multiply(v1, v2, mod1);
        auto mod_2_res = mod_multiply(v1, v2, mod2);
        assert(mod_1_res.size() == len && mod_2_res.size() == len);

        for (int i = 0; i < len; i++){
            mod_1_res[i] = mod_1_res[i] + (mod_2_res[i] - mod_1_res[i] + mod2) * mod1 % mod2 * mod1;
        }

        return mod_1_res;
    }

    /***
     * Computes the circular convolution of v1 and v2
     * The vectors must be of the same size, if not normalize and pad with zeros
     * Example to demonstrate convolution for n = 5, let A=v1, B=v2 and the result be C
     *
     * C0 = A0B0 + A1B4 + A2B3 + A3B2 + A4B1
     * C1 = A0B1 + A1B0 + A2B4 + A3B3 + A4B2
     * ...
     * ...
     * C4 = A0B4 + A1B3 + A2B2 + A3B1 + A4B0
     *
     * Note: If linear convolution is required (i.e, no wrapping around), pad with zeros accordingly
     *
     * Also read the notes on multiply()
     *
    ***/
    vector<long long> convolution(const vector <long long>& v1, const vector <long long>& v2){
        auto p = build_convolution(v1, v2);
        return multiply(p.first, p.second);
    }

    /***
     * Same as convolution(v1, v2), only values are calculated modulo mod
     *
    ***/
    vector<long long> mod_convolution(const vector <long long>& v1, const vector <long long>& v2, int mod){
        auto p = build_convolution(v1, v2);
        return mod_multiply(p.first, p.second, mod);
    }

    /***
     * Same as convolution(v1, v2), but supports larger values
     *
     * Also read the notes on ll_multiply()
     *
    ***/
    vector<long long> ll_convolution(const vector <long long>& v1, const vector <long long>& v2){
        auto p = build_convolution(v1, v2);
        return ll_multiply(p.first, p.second);
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
    vector<long long> hamming_distance(const char* str, const char* pattern){
        assert(is_binary_string(str) && is_binary_string(pattern));

        int n = strlen(str), m = strlen(pattern);
        vector<long long> v1(n, 0), v2(m, 0);

        for (int i = 0; i < n; i++) v1[i] = str[i] == '1' ? 1 : -1;
        for (int i = 0, j = m - 1; j >= 0; i++, j--) v2[i] = pattern[j] == '1' ? 1 : -1;

        vector<long long> res;
        auto v = multiply(v1, v2);
        for (int i = 0; (i + m) <= n; i++){
            res.push_back(m - ((A[i + m - 1] + m) >> 1));
        }
        return res;
    }

    /***
     * And convolution vector with every substring of length |pattern| in str (sharing only common 1 bits)
     *
     * str and pattern consists of only '1' and '0'
     *
     * For example,
     *     str: "0110110"
     *     pattern: "110"
     *     returns: [1, 2, 1, 1, 2]
     *
    ***/
    vector<long long> and_convolution(const char* str, const char* pattern){
        assert(is_binary_string(str) && is_binary_string(pattern));

        int n = strlen(str), m = strlen(pattern);
        vector<long long> v1(n, 0), v2(m, 0);

        for (int i = 0; i < n; i++) v1[i] = str[n - i - 1] - 48;
        for (int i = 0; i < m; i++) v2[i] = pattern[i] - 48;

        vector<long long> res;
        auto v = multiply(v1, v2);
        for (int i = n; i >= m; i--) res.push_back(v[i - 1]);
        return res;
    }
}

int main(){
    vector <long long> v1, v2, expected_result;

    v1 = {5, 1, 2, 6, 9, 8};
    v2 = {3, 9, 0, 2};
    assert(fft::multiply(v1, v2) == vector<long long>({15, 48, 15, 46, 83, 109, 84, 18, 16}));

    int mod = 14;
    assert(fft::mod_multiply(v1, v2, mod) == vector<long long>({1, 6, 1, 4, 13, 11, 0, 4, 2}));

    for (auto && x: v1) x = (1 << 30) - x;
    for (auto && x: v2) x = (1 << 30) - x;
    expected_result = {1152921496016912399, 55842988386341168, 1208764490845704463, 111685985362616878, 111685981067649619, 111685973551456877, 1208764475813318996, 55842987312599314, 1152921493869428752};

    assert(fft::multiply(v1, v2) != expected_result);  /// should fail because of precision, even with long double
    assert(fft::ll_multiply(v1, v2) == expected_result);

    v1 = {1, 2, 3, 4};
    v2 = {1, 0, 0, 2};
    assert(fft::convolution(v1, v2) == vector<long long>({1, 2, 3, 6, 5, 8, 11, 6, 4, 6, 8, 0, 0, 0, 0}));

    mod = 2;
    assert(fft::mod_convolution(v1, v2, mod) == vector<long long>({1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}));

    for (auto && x: v1) x = (1 << 30) - x;
    for (auto && x: v2) x = (1 << 30) - x;
    expected_result = {1152921502459363329, 55843003418726658, 1208764504804348163, 111686001468744198, 111686001468744197, 111686001468744200, 111686001468744203, 111686001468744198, 1208764500509380868, 55842998050017542, 1152921498164396040, 0, 0, 0, 0};

    assert(fft::convolution(v1, v2) != expected_result);  /// should fail because of precision, even with long double
    assert(fft::ll_convolution(v1, v2) == expected_result);

    expected_result = {3, 3, 1, 1, 4, 1, 2};
    assert(fft::hamming_distance("1000100101", "0110") == expected_result);

    expected_result = {1, 2, 1, 1, 2};
    assert(fft::and_convolution("0110110", "110") == expected_result);

    return 0;
}
