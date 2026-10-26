/***
 *
 * Digits of Factorial
 * Computes the number of digits of n! in base b
 *
 * Complexity: O(1)
 *
 * digits_of_factorial(n, b = 10): the digit count of n! in base b, for n >= 0 and b >= 2
 *
 * Notes:
 *   - Exact up to about n = 1e11, long double rounding gives wrong counts beyond
 *   - Below that too, b near a root of n! (b^k within long double precision of n!) can be off by one when n > 20
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long digits_of_factorial(long long n, long long b=10){
    assert(n >= 0 && b > 1);

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
    /// 10! = 3628800 has 7 digits, 100! has 158 and (10^9)! has 8565705523
    assert(digits_of_factorial(0) == 1);
    assert(digits_of_factorial(10) == 7);
    assert(digits_of_factorial(100) == 158);
    assert(digits_of_factorial(1000000000) == 8565705523);

    /// 5! = 120 is 1111000 in base 2, and 158 decimal digits make 79 in base 100
    assert(digits_of_factorial(5, 2) == 7);
    assert(digits_of_factorial(100, 100) == 79);

    return 0;
}
