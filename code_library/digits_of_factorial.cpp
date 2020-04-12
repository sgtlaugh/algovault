#include <stdio.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-11

const long double TWICE_PI = 4.0 * acos(0.0);
const long double HALF_LOG_TWICE_PI = 0.5 * log(TWICE_PI);

long long digits_of_factorial(long long n, long long base=10){
    assert(base > 1);
    if (n == 0 || n == 1) return 1;
    long double x = ((n + 0.5) * log(n) - n + HALF_LOG_TWICE_PI) / log(base);
    return ceil(x - EPS);
}

int main(){
    printf("%lld\n", digits_of_factorial(0));                        /// 1
    printf("%lld\n", digits_of_factorial(1));                        /// 1
    printf("%lld\n", digits_of_factorial(10));                       /// 7
    printf("%lld\n", digits_of_factorial(100));                      /// 158
    printf("%lld\n", digits_of_factorial(1000));                     /// 2568
    printf("%lld\n", digits_of_factorial(1000000));                  /// 5565709
    printf("%lld\n", digits_of_factorial(1000000000));               /// 8565705523

    printf("%lld\n", digits_of_factorial(100, 2));                   /// 525
    printf("%lld\n", digits_of_factorial(100, 100));                 /// 79
    printf("%lld\n", digits_of_factorial(2000000000, 666666667));    /// 2009706986
    printf("%lld\n", digits_of_factorial(1000000000, 2000000000));   /// 920941609

    return 0;
}
