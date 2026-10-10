/***
 *
 * Fast Walsh Hadamard Transformation to calculate convolution between two vectors
 * Convolution type can be either xor/or/and, exact over integers or modulo mod
 * Complexity for all convolutions: O(n log n)
 *
 * Convolutions of two vectors A, B of length N can be defined as below:
 * C = [0] * n
 * for i in range(0, n):
 *     for j in range(0, n):
 *          C[i operator j] += A[i] * B[j] // where operator is one of {^, |, &}
 *
 * Notes:
 *   - A and B must be of the same length n
 *   - n must be a power of 2
 *   - Exact mode (mod = 0): no overflow while n^2 * max|A[i]| * max|B[j]| < 2^62
 *   - Modular mode: 1 <= mod < 2^31, any input values (negative or >= mod are reduced), results in [0, mod)
 *   - xor needs an odd mod, the inverse halves each level by multiplying with inv(2); or/and take any mod
 *   - Calling walsh_transform / inverse_walsh_transform directly with a mod requires values already in [0, mod)
 *   - No length cap; use one instance per thread (it owns a scratch buffer)
 *
 * Usage:
 *   WalshHadamard fwht;
 *   vector<long long> C = fwht.xor_convolution(A, B);             // also or_convolution, and_convolution
 *   vector<long long> D = fwht.xor_convolution(A, B, 998244353);  // modular mode
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct WalshHadamard{
    static constexpr int OR = 0;
    static constexpr int AND = 1;
    static constexpr int XOR = 2;

    vector<long long> buffer;  /// reused across calls: allocating a fresh 8 MB vector per call at n = 2^20 costs ~10% in page faults

    static void walsh_transform(long long* ar, int n, int conv_type, long long mod = 0){
        dispatch<false>(ar, n, conv_type, mod);
    }

    static void inverse_walsh_transform(long long* ar, int n, int conv_type, long long mod = 0){
        dispatch<true>(ar, n, conv_type, mod);
    }

    vector<long long> convolution(const vector<long long>& A, const vector<long long>& B, int conv_type, long long mod = 0){
        int n = A.size();
        assert(A.size() == B.size() && __builtin_popcount(n) == 1);
        vector<long long> res = A;
        buffer.assign(B.begin(), B.end());
        if (mod){
            for (auto& x : res) if ((x %= mod) < 0) x += mod;
            for (auto& x : buffer) if ((x %= mod) < 0) x += mod;
        }

        walsh_transform(res.data(), n, conv_type, mod);
        walsh_transform(buffer.data(), n, conv_type, mod);
        if (mod) for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i] % mod;
        else for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i];
        inverse_walsh_transform(res.data(), n, conv_type, mod);

        return res;
    }

    vector<long long> or_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, OR, mod);
    }

    vector<long long> and_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, AND, mod);
    }

    vector<long long> xor_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, XOR, mod);
    }

private:
    /// modular operands are in [0, mod); (s >> 63) & mod adds mod back exactly when s went negative
    template<bool modular>
    static long long add(long long x, long long y, long long mod){
        long long s = x + y;
        if constexpr (modular) s -= mod, s += (s >> 63) & mod;
        return s;
    }

    /// conv_type and mod != 0 become template arguments so the butterfly loops carry no per-element branches
    template<bool inverse>
    static void dispatch(long long* ar, int n, int conv_type, long long mod){
        assert(mod >= 0 && mod < (1LL << 31) && (conv_type != XOR || !mod || mod % 2 == 1));
        if (conv_type == OR) mod ? transform<OR, inverse, true>(ar, n, mod) : transform<OR, inverse, false>(ar, n, mod);
        if (conv_type == AND) mod ? transform<AND, inverse, true>(ar, n, mod) : transform<AND, inverse, false>(ar, n, mod);
        if (conv_type == XOR) mod ? transform<XOR, inverse, true>(ar, n, mod) : transform<XOR, inverse, false>(ar, n, mod);
    }

    /// exact mode: v is even, so >> 1 divides; modular mode: mod is odd, so adding it to an odd v makes v + mod even
    template<bool modular>
    static long long half(long long v, long long mod){
        if constexpr (modular) v += (v & 1) * mod;
        return v >> 1;
    }

    template<bool modular>
    static long long sub(long long x, long long y, long long mod){
        long long d = x - y;
        if constexpr (modular) d += (d >> 63) & mod;
        return d;
    }

    template<int conv_type, bool inverse, bool modular>
    static void transform(long long* ar, int n, long long mod){
        if (n < 2) return;

        int m = n >> 1;
        transform<conv_type, inverse, modular>(ar, m, mod);
        transform<conv_type, inverse, modular>(ar + m, m, mod);

        for (int i = 0; i < m; i++){
            long long x = ar[i], y = ar[i + m];
            if constexpr (conv_type == OR) ar[i + m] = inverse ? sub<modular>(y, x, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == AND) ar[i] = inverse ? sub<modular>(x, y, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == XOR){
                ar[i] = add<modular>(x, y, mod), ar[i + m] = sub<modular>(x, y, mod);
                if constexpr (inverse) ar[i] = half<modular>(ar[i], mod), ar[i + m] = half<modular>(ar[i + m], mod);
            }
        }
    }
};

int main(){
    WalshHadamard fwht;
    const vector<long long> A = {0, 1, 3, 5};
    const vector<long long> B = {1, 1, 2, 1};

    assert(fwht.or_convolution(A, B) == vector<long long>({0, 2, 9, 34}));
    assert(fwht.and_convolution(A, B) == vector<long long>({14, 7, 19, 5}));
    assert(fwht.xor_convolution(A, B) == vector<long long>({12, 14, 9, 10}));

    assert(fwht.or_convolution(A, B, 7) == vector<long long>({0, 2, 2, 6}));
    assert(fwht.and_convolution(A, B, 4) == vector<long long>({2, 3, 3, 1}));
    assert(fwht.xor_convolution(A, B, 7) == vector<long long>({5, 0, 2, 3}));

    assert(fwht.or_convolution({-7}, {6}) == vector<long long>({-42}));
    assert(fwht.and_convolution({-7}, {6}) == vector<long long>({-42}));
    assert(fwht.xor_convolution({-7}, {6}) == vector<long long>({-42}));
    assert(fwht.xor_convolution({-7}, {6}, 5) == vector<long long>({3}));
    assert(fwht.xor_convolution({1000000000000000000LL, 3}, {1, 1}, 1) == vector<long long>({0, 0}));

    const int n = 1 << 21;
    vector<long long> big(n), delta_low(n), delta_high(n), expected_and(n), expected_xor(n);
    for (int i = 0; i < n; i++) big[i] = i % 7 - 3;
    delta_low[0] = 1, delta_high[n - 1] = 2;
    for (int i = 0; i < n; i++) expected_and[i] = 2 * big[i], expected_xor[i] = 2 * big[i ^ (n - 1)];

    assert(fwht.or_convolution(big, delta_low) == big);
    assert(fwht.and_convolution(big, delta_high) == expected_and);
    assert(fwht.xor_convolution(big, delta_high) == expected_xor);

    const long long mod = 998244353;
    const int m = 1 << 20;
    vector<long long> all_max(m, mod - 1);
    assert(fwht.xor_convolution(all_max, all_max, mod) == vector<long long>(m, m));

    return 0;
}
