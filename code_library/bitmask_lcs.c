#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define MAX 50010


bool flag[MAX];
char A[MAX], B[MAX], S[2][MAX];


/***
 * LCS with bitmasking
 * Complexity: O (N * M / 64)
***/

int lcs(char* A, char* B){
    int i, j, k, v, n, m, res = 0;
    unsigned long long x, t, q, y, mask[128] = {0};

    memset(flag, 0, sizeof(flag));
    for (n = 0; A[n]; n++) S[0][n] = A[n];
    for (m = 0; B[m]; m++) S[1][m] = B[m];

    for (i = 0; (i * 64) < m; i++){
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
    srand(666);

    int i, n, m;
    n = MAX - 10, m = MAX - 10;
    for (i = 0; i < n; i++) A[i] = (rand() % 26) + 'a';
    for (i = 0; i < m; i++) B[i] = (rand() % 26) + 'a';
    A[n] = B[m] = 0;

    clock_t start = clock();
    printf("%d\n", lcs(A, B));  /// 16245

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// took 0.3 s locally
    return 0;
}
