#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#define MAXK 1010
#define MOD 1000000007

/***
 *
 * Faulhaber's formula to calculate the sum of k-th powers
 *
 * Example:
 *   - faulhaber(5, 2) = 1^2 + 2^2 + 3^2 + 4^2 + 5^2 = 55
 *
 * Complexity:
 *   - O(MAXK^2) pre-processing once to compute stirling numbers
 *   - O(k) for each call to compute sum of powers
 *
***/

int S[MAXK][MAXK], inv[MAXK];

int expo(long long x, int n){
    x %= MOD;
    long long res = 1;

    while (n){
        if (n & 1) res = res * x % MOD;
        x = x * x % MOD;
        n >>= 1;
    }

    return res % MOD;
}

void generate(){
    int i, j;
    for (i = 0; i < MAXK; i++) inv[i] = expo(i, MOD - 2);

    S[0][0] = 1;
    for (i = 1; i < MAXK; i++){
        S[i][0] = 0;
        for (j = 1; j <= i; j++){
            S[i][j] = ((long long)S[i - 1][j] * j + S[i - 1][j - 1]) % MOD;
        }
    }
}

int faulhaber(long long n, int k){
    n %= MOD;
    if (!k) return n;

    int j;
    long long res = 0, p = 1;
    for (j = 0; j <= k; j++){
        p = p * (n + 1 - j) % MOD;
        res = (res + (S[k][j] * p % MOD) * inv[j + 1]) % MOD;
    }

    return res % MOD;
}

int main(){
    generate();

    assert(faulhaber(1, 1) == 1);
    assert(faulhaber(9, 0) == 9);
    assert(faulhaber(5, 2) == 55);
    assert(faulhaber(1000000000, 512) == 244343556);

    return 0;
}
