/***
 *
 * Computes the number of digits of n! in base b
 *
 * digits_of_factorial(n, b = 10): the digit count of n! in base b, for n >= 0 and b >= 2
 * Exact for n <= 20, uses lgammal beyond that
 * Although fairly accurate, there is no guarantee it'll compute the exact answer for larger n
 *
 * Complexity: O(1)
 *
***/

#include <math.h>
#include <assert.h>

long long digits_of_factorial(long long n, long long b=10){
    assert(b > 1);

    /// n! fits in 64 bits, an exact count handles b = n!, the only exact power and one a logarithm can round either way
    if (n <= 20){
        unsigned long long f = 1;
        long long res = 0;
        for (int i = 2; i <= n; i++) f *= i;
        for (; f; f /= b) res++;
        return res;
    }

    return floorl(lgammal(n + 1.0L) / logl(b)) + 1;
}

#include <stdio.h>

int main(){
    int t;
    scanf("%d", &t);
    while (t--){
        int n;
        scanf("%d", &n);
        printf("%lld\n", digits_of_factorial(n));
    }
    return 0;
}
