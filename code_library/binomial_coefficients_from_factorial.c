#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define MAX 1000010
#define MOD 1000000007

int factorial[MAX];

int expo(int x, int n){
    long long res = 1;

    while (n){
        if (n & 1) res = res * x % MOD;
        x = (long long)x * x % MOD;
        n >>= 1;
    }

    return res % MOD;
}

/// binomial(n, k) modulo MOD
int binomial(int n, int k){
    int v = (long long)factorial[k] * factorial[n - k] % MOD;
    return (long long)factorial[n] * expo(v, MOD - 2) % MOD;
}

void generate(){
    int i;
    for (i = 1, factorial[0] = 1; i < MAX; i++){
        factorial[i] = (long long)factorial[i - 1] * i % MOD;
    }
}

int main(){
    generate();

    printf("%d\n", binomial(1, 0));           /// 1
    printf("%d\n", binomial(1, 1));           /// 1
    printf("%d\n", binomial(10, 5));          /// 252
    printf("%d\n", binomial(1000000, 99966)); /// 978206100

    return 0;
}
