/***
 *
 * Computes the number of digits of n! in base b
 * Exact for n <= 20, uses lgammal beyond that
 * Although fairly accurate, there is no guarantee it'll compute the exact answer for larger n
 *
***/

#include <stdio.h>
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

int main(){
    assert(digits_of_factorial(0) == 1);
    assert(digits_of_factorial(1) == 1);
    assert(digits_of_factorial(2) == 1);
    assert(digits_of_factorial(3) == 1);
    assert(digits_of_factorial(4) == 2);
    assert(digits_of_factorial(5) == 3);
    assert(digits_of_factorial(10) == 7);
    assert(digits_of_factorial(100) == 158);
    assert(digits_of_factorial(1000) == 2568);
    assert(digits_of_factorial(1000000) == 5565709);
    assert(digits_of_factorial(1000000000) == 8565705523);

    assert(digits_of_factorial(100, 2) == 525);
    assert(digits_of_factorial(100, 100) == 79);
    assert(digits_of_factorial(2000000000, 666666667) == 2009706986);
    assert(digits_of_factorial(1000000000, 2000000000) == 920941609);

    return 0;
}
