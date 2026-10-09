#include "../common.h"

#define main library_main
#include "../../code_library/misc/divisors.cpp"
#undef main

/// Divisors by testing every candidate up to sqrt(n)
vector<long long> brute(long long n){
    vector<long long> small, large;
    for (long long d = 1; d * d <= n; d++){
        if (n % d) continue;
        small.push_back(d);
        if (d * d != n) large.push_back(n / d);
    }
    small.insert(small.end(), large.rbegin(), large.rend());
    return small;
}

int main(){
    for (long long n = 1; n <= 20000; n++) assert(divisors(n) == brute(n));
    for (long long it = 0; it < stress::scaled(300); it++){
        long long n = it % 3 == 0 ? stress::rand_int(1, 100000000000000LL) : stress::rand_int(1, 1000000000000LL);
        assert(divisors(n) == brute(n));
    }
    /// Highly composite numbers and a prime square at the 1e14 bound
    for (long long n : {963761198400LL, 97821761637600LL, 9999991LL * 9999991LL, 100000000000000LL}) assert(divisors(n) == brute(n));
    return 0;
}
