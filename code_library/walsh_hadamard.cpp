#include <bits/stdc++.h>

#define MAX 1048576

using namespace std;

namespace fwht{
    const int OR =  0;
    const int AND = 1;
    const int XOR = 2;

    long long P1[MAX], P2[MAX];

    void walsh_transform(long long* ar, int n, int conv_type){
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

    void inverse_walsh_transform(long long* ar, int n, int conv_type){
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
        for (int i = 0; i < n; i++) P1[i] = A[i];
        for (int i = 0; i < n; i++) P2[i] = B[i];

        walsh_transform(P1, n, conv_type);
        walsh_transform(P2, n, conv_type);
        for (int i = 0; i < n; i++) P1[i] = P1[i] * P2[i];
        inverse_walsh_transform(P1, n, conv_type);
        return vector<long long> (P1, P1 + n);
    }


    /***
     * C = [0] * n
     * for i in range(0, n):
     *     for j in range(0, n):
     *          C[i | j] += A[i] * B[j]
     *
     * return C as a vector in n log n, n must be a power of 2
     *
    ***/
    vector <long long> or_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, OR);
    }

    /***
     * C = [0] * n
     * for i in range(0, n):
     *     for j in range(0, n):
     *          C[i & j] += A[i] * B[j]
     *
     * return C as a vector in n log n, n must be a power of 2
     *
    ***/
    vector <long long> and_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, AND);
    }

    /***
     * C = [0] * n
     * for i in range(0, n):
     *     for j in range(0, n):
     *          C[i ^ j] += A[i] * B[j]
     *
     * return C as a vector in n log n, n must be a power of 2
     *
    ***/
    vector <long long> xor_convolution(const vector <long long>& A, const vector <long long>& B){
        return convolution(A, B, XOR);
    }
}

int main(){
    const vector <long long> A = {0, 1, 3, 5};
    const vector <long long> B = {1, 1, 2, 1};
    
    assert(fwht::or_convolution(A, B) == vector<long long>({0, 2, 9, 34}));
    assert(fwht::and_convolution(A, B) == vector<long long>({14, 7, 19, 5}));
    assert(fwht::xor_convolution(A, B) == vector<long long>({12, 14, 9, 10}));

    return 0;
}
