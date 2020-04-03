#include <stdio.h>

#define MAX 3010
#define MOD 1000000007

/***
 *
 * https://en.wikipedia.org/wiki/Bernoulli_number
 *
 * computes all the bernoulli numbers from 1 to N in O(N^2)
 * bernoulli numbers are sequence of rational numbers of the form p/q
 * since the numbers can be huge, they are calculated modulo MOD
 * bernoulli[i] = p * inv(q, MOD) % MOD
 *
***/


int S[MAX][MAX], inv[MAX], factorial[MAX], bernoulli[MAX];

int expo(long long x, int n){
    long long res = 1;

    while (n){
        if (n & 1) res = (res * x) % MOD;
        x = (x * x) % MOD;
        n >>= 1;
    }

    return (res % MOD);
}

void generate(){
    int i, j;
    long long x, y, z, lim = (long long)MOD * MOD;

    for (i = 1, factorial[0] = 1; i < MAX; i++) factorial[i] = ((long long) factorial[i - 1] * i) % MOD;
    for (i = 0; i < MAX; i++) inv[i] = expo(i, MOD - 2);
    for (i = 1, S[0][0] = 1; i < MAX; i++){
        for (j = 1, S[i][0] = 0; j <= i; j++) S[i][j] = ( ((long long)S[i - 1][j] * j) + S[i - 1][j - 1]) % MOD;
    }

    bernoulli[0] = 1;
    for (i = 1; (i + 1) < MAX; i++){
        if ((i & 1) && i > 1) bernoulli[i] = 0;
        else{
            for (j = 0, x = 0, y = 0; j <= i; j++){
                z = ((long long) factorial[j] * inv[j + 1]) % MOD;
                z = (z * S[i][j]) % MOD;
                if (j & 1) y += z;
                else x += z;
            }
            bernoulli[i] = (lim + x - y) % MOD;
        }
    }
}

int main(){
    generate();
    printf("%d\n", bernoulli[10]);  /// bernoulli[10] = 5 / 66 = 5 * 469696973 % 1000000007 = 348484851
    return 0;
}
