#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define MAX 1000010

int n, ar[MAX];

/// https://en.wikipedia.org/wiki/Comb_sort
void combsort(int n, int* ar){
    if (n <= 1) return;
    int i, x, g = n, flag = 0;
    while ((g != 1) || flag){
        flag = 0;
        if (g != 1) g *= 0.77425;

        for (i = 0; (i + g) < n; i++){
            if (ar[i] > ar[i + g]){
                flag = 1;
                x = ar[i], ar[i] = ar[i + g], ar[i + g] = x;
            }
        }
    }
}

int main(){
    int i;
    n = MAX - 10;
    for (i = 0; i < n; i++) ar[i] = (rand() << 15 ^ rand()) % 1000000007;

    clock_t start = clock();
    combsort(n, ar);
    fprintf(stderr, "Time taken = %0.5f s\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 0.232 seconds locally

    for (i = 0; (i + 1) < n; i++) assert(ar[i] <= ar[i + 1]);
    return 0;
}
