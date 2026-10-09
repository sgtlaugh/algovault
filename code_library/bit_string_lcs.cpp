/***
 *
 * Bit-String LCS (http://users.monash.edu/~lloyd/tildeStrings/Alignment/86.IPL.html)
 * Calculates the longest common subsequence of two strings with bitset
 *
 * Complexity: O(N * M / 64)
 *
***/

#include <bits/stdc++.h>

#define MAX 100000

using namespace std;

/// Hyyro's row update V = (V + U) | (V - U) with U = V & match, zeros of V count the LCS
/// U is a submask of V so V - U never borrows, only the addition carries into the next 64-bit block
int lcs(const char* A, const char* B){
    int i, j, k, n = strlen(A), m = strlen(B), res = 0;
    unsigned long long u, v, sum, mask[256];
    vector<char> carry(n, 0);

    for (i = 0; i * 64 < m; i++){
        memset(mask, 0, sizeof(mask));
        for (k = 0; k < 64 && i * 64 + k < m; k++){
            mask[(unsigned char)B[i * 64 + k]] |= (1ULL << k);
        }

        /// Bits past m stay set, since mask is 0 there and V - U keeps them
        for (j = 0, v = ~0ULL; j < n; j++){
            u = v & mask[(unsigned char)A[j]];
            bool c = __builtin_uaddll_overflow(v, u, &sum);
            c |= __builtin_uaddll_overflow(sum, carry[j], &sum);
            carry[j] = c;
            v = sum | (v - u);
        }

        res += __builtin_popcountll(~v);
    }
    return res;
}

int main(){
    char A[MAX], B[MAX];

    int n = MAX - 10, m = MAX - 10;
    for (int i = 0; i < n; i++) A[i] = ((long long)i * i % 26) + 'a';
    for (int i = 0; i < m; i++) B[i] = ((long long)i * i * i % 26) + 'a';
    A[n] = B[m] = 0;

    clock_t start = clock();
    assert(lcs(A, B) == 23075);

    fprintf(stderr, "Time taken = %0.3f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// 0.715 s locally
    return 0;
}
