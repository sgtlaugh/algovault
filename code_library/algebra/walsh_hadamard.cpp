/***
 *
 * Fast Walsh Hadamard Transformation to calculate convolution between two vectors
 * Convolution type can be either xor/or/and
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
 *   - No overflow while n^2 * max|A[i]| * max|B[j]| < 2^62
 *   - No length cap; use one instance per thread (it owns a scratch buffer)
 *
 * Usage:
 *   WalshHadamard fwht;
 *   vector <long long> C = fwht.xor_convolution(A, B);  // also or_convolution, and_convolution
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct WalshHadamard{
    static constexpr int OR = 0;
    static constexpr int AND = 1;
    static constexpr int XOR = 2;

    vector <long long> buffer;  /// reused across calls: allocating a fresh 8 MB vector per call at n = 2^20 costs ~10% in page faults

    static void walsh_transform(long long* ar, int n, int conv_type){
        if (!n) return;

        int m = n >> 1;
        walsh_transform(ar, m, conv_type);
        walsh_transform(ar + m, m, conv_type);

        /// Add modulo operations below if required
        for (int i = 0; i < m; i++){
            long long x = ar[i], y = ar[i + m];
            if (conv_type == OR) ar[i] = x, ar[i + m] = x + y;
            if (conv_type == AND) ar[i] = x + y, ar[i + m] = y;
            if (conv_type == XOR) ar[i] = x + y, ar[i + m] = x - y;
        }
    }

    static void inverse_walsh_transform(long long* ar, int n, int conv_type){
        if (!n) return;

        int m = n >> 1;
        inverse_walsh_transform(ar, m, conv_type);
        inverse_walsh_transform(ar + m, m, conv_type);

        /// Add modulo operations below if required, inverse modulo may be requried for XOR
        for (int i = 0; i < m; i++){
            long long x = ar[i], y = ar[i + m];
            if (conv_type == OR) ar[i] = x, ar[i + m] = y - x;
            if (conv_type == AND) ar[i] = x - y, ar[i + m] = y;
            if (conv_type == XOR) ar[i] = (x + y) >> 1, ar[i + m] = (x - y) >> 1;
        }
    }

    vector <long long> convolution(const vector <long long>& A, const vector <long long>& B, int conv_type){
        int n = A.size();
        assert(A.size() == B.size() && __builtin_popcount(n) == 1);
        vector <long long> res = A;
        buffer.assign(B.begin(), B.end());

        walsh_transform(res.data(), n, conv_type);
        walsh_transform(buffer.data(), n, conv_type);
        for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i];
        inverse_walsh_transform(res.data(), n, conv_type);

        return res;
    }

    vector <long long> or_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, OR);
    }

    vector <long long> and_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, AND);
    }

    vector <long long> xor_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, XOR);
    }
};

int main(){
    WalshHadamard fwht;
    const vector <long long> A = {0, 1, 3, 5};
    const vector <long long> B = {1, 1, 2, 1};

    assert(fwht.or_convolution(A, B) == vector<long long>({0, 2, 9, 34}));
    assert(fwht.and_convolution(A, B) == vector<long long>({14, 7, 19, 5}));
    assert(fwht.xor_convolution(A, B) == vector<long long>({12, 14, 9, 10}));

    assert(fwht.or_convolution({-7}, {6}) == vector<long long>({-42}));
    assert(fwht.and_convolution({-7}, {6}) == vector<long long>({-42}));
    assert(fwht.xor_convolution({-7}, {6}) == vector<long long>({-42}));

    const int n = 1 << 21;
    vector <long long> big(n), delta_low(n), delta_high(n), expected_and(n), expected_xor(n);
    for (int i = 0; i < n; i++) big[i] = i % 7 - 3;
    delta_low[0] = 1, delta_high[n - 1] = 2;
    for (int i = 0; i < n; i++) expected_and[i] = 2 * big[i], expected_xor[i] = 2 * big[i ^ (n - 1)];

    assert(fwht.or_convolution(big, delta_low) == big);
    assert(fwht.and_convolution(big, delta_high) == expected_and);
    assert(fwht.xor_convolution(big, delta_high) == expected_xor);

    return 0;
}
