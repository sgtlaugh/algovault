/***
 *
 * Bit-String LCS (http://users.monash.edu/~lloyd/tildeStrings/Alignment/86.IPL.html)
 * Calculates the longest common subsequence of two strings with bitset
 *
 * Complexity: O(N * M / 64)
 *
***/

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#define MAX 50010

bool flag[MAX];
char A[MAX], B[MAX], S[2][MAX];

int lcs(char* A, char* B){
    int i, j, k, v, n, m, res = 0;
    unsigned long long x, t, q, y, mask[128] = {0};

    memset(flag, 0, sizeof(flag));
    for (n = 0; A[n]; n++) S[0][n] = A[n];
    for (m = 0; B[m]; m++) S[1][m] = B[m];

    for (i = 0; i * 64 < m; i++){
        memset(mask, 0, sizeof(mask));
        for (k = 0; k < 64 && (i * 64 + k) < m; k++){
            mask[(int)S[1][i * 64 + k]] |= (1ULL << k);
        }

        for (j = 0, x = 0; j < n; j++){
            t = mask[(int)S[0][j]] & ~x;
            x |= t, v = flag[j];
            q = x - (t << 1) - v, y = (q & ~x) | t;
            flag[j] = y >> 63, x &= ~(y << 1);
            if (v) x &= ~1ULL;
        }
        res += __builtin_popcountll(x);
    }
    return res;
}

int main(){
    int n = MAX - 10, m = MAX - 10;
    for (int i = 0; i < n; i++) A[i] = ((long long)i * i % 26) + 'a';
    for (int i = 0; i < m; i++) B[i] = ((long long)i * i * i % 26) + 'a';
    A[n] = B[m] = 0;

    clock_t start = clock();
    assert(lcs(A, B) == 11540);

    fprintf(stderr, "Time taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// 0.234000 s locally
    return 0;
}
