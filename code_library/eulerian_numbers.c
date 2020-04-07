#include <stdio.h>

#define MAX 5010
#define MOD 1000000007

/***

Eulerian number A(n, k) is the number of permutations of 1 to n in which exactly k elements are greater than their previous element
Eulerian triangle for n = 1 to 7 and k = 0 to n - 1 below


1
1 1
1 4 1
1 11 11 1
1 26 66 26 1
1 57 302 302 57 1
1 120 1191 2416 1191 120 1

***/

int dp[MAX][MAX];

void generate(){
    int n, k;
    for (n = 0; n < MAX; n++){
        for (k = 1, dp[n][0] = 1; k <= n; k++){
            dp[n][k] = ((long long)dp[n - 1][k] * (k + 1) + (long long)dp[n - 1][k - 1] * (n - k)) % MOD;
        }
    }
}

int main(){
    int n, k;
    generate();

    for (n = 1; n <= 7; n++){
        for (k = 0; k < n; k++){
            printf("%d ", dp[n][k]);
        }
        puts("");
    }

    return 0;
}
