/***

 * USACO 2017 Jan Platinum 3 - SUBSEQUENCE REVERSAL (https://goo.gl/y39FAY)

 * Given an array, it is possible to reverse any subsequence at most once
 * Find the length of the longest non-decreasing sequence after using at most one reversal

 * IMPORTANT: p = exp((double)-(next - cur) / temperature) when minimizing function

***/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX 52
#define clr(ar) memset(ar, 0, sizeof(ar))


const double multiplier = 0.997; /// Usually good to have between 0.95 - 0.999


int n, m, dp[MAX], ar[MAX], adj[MAX], bits[MAX], last_bits[MAX];

double rand_real(){
    return (double)rand() / RAND_MAX;
}

int solve(){
    int i, j, x, res = 0;

    for (i = n - 1; i >= 0; i--){
        dp[i] = 1;
        for (j = i + 1; j < n; j++){
            if (ar[i] <= ar[j]){
                x = 1 + dp[j];
                if (x > dp[i]) dp[i] = x;
            }
        }
        if (dp[i] > res) res = dp[i];
    }
    return res;
}

int main(){
    srand(17);
    int t, i, j, k, x, a, b, lim, len, cur, res, next;

    while (scanf("%d", &n) != EOF){
        for (i = 0; i < n; i++) scanf("%d", &ar[i]);

        clr(bits), clr(last_bits);
        cur = res = solve(), lim = 180000; /// lim = total iterations, the more the better

        double temperature = 1e9; /// initial temperature
        for (t = 0; t < lim; t++){
            /// m = number of indexes selected for reversal
            m = ((lim - t - 1) * n / lim) + 1; /// select more during first and gradually decrease eventually for stability
            for (i = 0; i < n; i++) last_bits[i] = bits[i];

            while (m--) bits[rand() % n] ^= 1;
            for (i = 0, len = 0; i < n; i++){
                if (bits[i]) adj[len++] = i;
            }

            for (k = 0; k < 2; k++){
                for (i = 0, j = len - 1; i < j; i++, j--){
                    a = adj[i], b = adj[j];
                    x = ar[a], ar[a] = ar[b], ar[b] = x;
                }
                if (k == 0) next = solve();
            }
            if (next > res) res = next; /// update global maximum

            if (next > cur) cur = next; /// if state is better than before update
            else{
                /// IMPORTANT: p = exp((double)-(next - cur) / temperature) when minimizing function
                double p = exp((double)(next - cur) / temperature);

                if (rand_real() <= p) cur = next; /// otherwise update with a probability depending on the temperature
                else{
                    for (i = 0; i < n; i++) bits[i] = last_bits[i];
                }
            }

            if (temperature > 1e-20) temperature *= multiplier; /// update temperature after each iteration
        }

        printf("%d\n", res);
    }
    return 0;
}
