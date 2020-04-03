#include <stdio.h>
#include <string.h>

#define MAX 1010
#define MOD 1000000007

/// calculate nCk and nPk with DP
/// nCk[n][k] = binomial(n, k)
/// nPk[n][k] = nCk[n][k] * k!
int nCk[MAX][MAX], nPk[MAX][MAX];

void generate(){
    int i, j;
    memset(nCk, 0, sizeof(nCk));

    for (i = 0; i < MAX; i++){
        for (j = 0; j <= i; j++){
            if (i == j || j == 0) nCk[i][j] = 1;
            else{
                nCk[i][j] = nCk[i - 1][j] + nCk[i - 1][j - 1];
                if (nCk[i][j] >= MOD) nCk[i][j] -= MOD;
            }
        }
    }

    memset(nPk, 0, sizeof(nPk));
    for (i = 0; i < MAX; i++){
        for (j = 0; j <= i; j++){
            if (!j) nPk[i][j] = 1;
            else nPk[i][j] = (long long)nPk[i][j - 1] * (i - j + 1) % MOD;
        }
    }
}

int main(){
    generate();

    printf("%d\n", nCk[10][3]);   /// 120
    printf("%d\n", nPk[10][3]);   /// 720

    return 0;
}
