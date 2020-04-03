#include <stdio.h>
#include <math.h>


const long double TWICE_PI = 4.0 * acos(0.0);
const long double LOG_10_E = log10(exp(1.0));
/*
/// Stirling's approximation to calculate the number of digits of n!
long long digits_of_factorial(long long n){
    if (n == 0 || n == 1) return 1;
    long double x = (0.5 * log(TWICE_PI * n) + (n * log(n)) - n) * LOG_10_E;
    return ceil(x);
}*/

long long digits_of_factorial(long long n, long long base){
    if (n == 0 || n == 1) return 1;
    long double x = (0.5 * log(TWICE_PI * n) + (n * log(n)) - n) * LOG_10_E;
    return ceil(x);
}

int main(){
    printf("%lld\n", digits_of_factorial(0));           /// 1
    printf("%lld\n", digits_of_factorial(1));           /// 1
    printf("%lld\n", digits_of_factorial(10));          /// 7
    printf("%lld\n", digits_of_factorial(100));         /// 158
    printf("%lld\n", digits_of_factorial(1000));        /// 2568
    printf("%lld\n", digits_of_factorial(1000000));     /// 5565709
    printf("%lld\n", digits_of_factorial(1000000000));  /// 8565705523

    return 0;
}
///18488885.0
