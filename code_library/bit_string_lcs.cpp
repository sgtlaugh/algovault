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

int lcs(const char* A, const char* B){
    bool flag[MAX];
    int i, j, k, v, n, m, res = 0;
    unsigned long long x, t, q, y, mask[128] = {0};

    memset(flag, 0, sizeof(flag));
    n = strlen(A), m = strlen(B);

    for (i = 0; i * 64 < m; i++){
        memset(mask, 0, sizeof(mask));
        for (k = 0; k < 64 && i * 64 + k < m; k++){
            mask[(int)B[i * 64 + k]] |= (1ULL << k);
        }

        for (j = 0, x = 0; j < n; j++){
            t = mask[(int)A[j]] & ~x;
            x |= t;
            v = flag[j];
            q = x - (t << 1) - v;
            y = (q & ~x) | t;
            flag[j] = y >> 63;
            x &= ~(y << 1);
            if (v) x &= ~1;
        }

        res += __builtin_popcountll(x);
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
