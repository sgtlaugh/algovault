/***
 * Calculates the mobius function for all values from 1 to MAX
 * 
 * mu[1] = 1, mu[n] = 0 if n has a squared prime factor
 * mu[n] = 1 if n is square-free with even number of prime factors
 * mu[n] = -1 if n is square-free with odd number of prime factors
 *
***/

#include <bits/stdc++.h>

using namespace std;

#define MAX 10000010
#define PRIME_CNT 666666  // Prime count from 1 to MAX

int len = 0, prime[PRIME_CNT];
char mu[MAX] = {0}, flag[MAX] = {0};

void generate_mobius_fast(){
    mu[1] = 1;
    int i, j, k;

    for (i = 2; i < MAX; i++){
        if (!flag[i]) mu[i] = -1, prime[len++] = i;
        for (j = 0; j < len && (k = i * prime[j]) < MAX; j++){
            flag[k] = 1;
            if (!(i % prime[j])){
                mu[k] = 0;
                break;
            }
            else mu[k] -= mu[i];
        }
    }
}

/// less code and elegant idea exploiting multiplicative functions, although slower
void generate_mobius(){
    mu[1] = 1;
    for (int i = 1; i < MAX; i++){
        for (int j = i + i; j < MAX; j += i){
            mu[j] -= mu[i];
        }
    }
}

int main(){
    generate_mobius_fast();

    int n;
    for (n = 0; n < 10; n++) printf("%d ", mu[n]);  /// 0 1 -1 -1 0 -1 1 -1 0 0
    puts("");

    return 0;
}
