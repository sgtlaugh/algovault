// LINK: -lquadmath
#include "common.h"
#include <quadmath.h>

#define main library_main
#include "../code_library/gamblers_ruin.c"
#undef main

/// Reference: solve P_i = p P_(i+1) + q P_(i-1), P_0 = 1, P_N = 0 with the Thomas algorithm (diagonally dominant, stable)
long double reference(int n1, int n2, long double p){
    if (n1 == 0) return 1;
    if (n2 == 0) return 0;
    int n = n1 + n2;
    long double q = 1 - p;
    vector<long double> c(n), d(n);  /// forward sweep over the unknowns P_1 .. P_(n-1), row i is -q P_(i-1) + P_i - p P_(i+1) = 0
    for (int i = 1; i < n; i++){
        long double denom = (i == 1) ? 1 : 1 + q * c[i - 1];
        c[i] = -p / denom;
        d[i] = (i == 1 ? q : q * d[i - 1]) / denom;  /// P_0 = 1 moves q into the first right hand side
    }
    long double x = 0;  /// back substitution from P_n = 0 down to P_n1
    for (int i = n - 1; i >= n1; i--) x = d[i] - c[i] * x;
    return x;
}

/// The same closed form in quad precision (113 bit mantissa), ~1e-30 accurate, far beyond the long double result
__float128 quad_reference(int n1, int n2, long double p_in){
    __float128 p = p_in, q = 1 - p, d = p - q, n = (__float128)n1 + n2;
    if (d == 0) return n2 / n;
    if (p < q){
        __float128 L = log1pq(d / q);
        return expm1q(n2 * L) / expm1q(n * L);
    }
    __float128 L = log1pq(-d / p);
    return expq(n1 * L) * expm1q(n2 * L) / expm1q(n * L);
}

int main(){
    auto close = [](long double a, long double b){ return fabsl(a - b) <= 1e-12L * max((long double)1, fabsl(b)); };

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n1 = stress::rand_int(0, 60), n2 = stress::rand_int(0, 60);
        if (n1 + n2 == 0) continue;
        long double p = stress::rand_int(0, 3) == 0 ? 0.5L + stress::rand_int(-1000, 1000) * 1e-15L : stress::rand_int(0, 1000000) / 1e6L;
        long double got = gamblers_ruin(n1, n2, p);
        assert(close(got, reference(n1, n2, p)));
        assert(0 <= got && got <= 1);
    }

    /// Large n against the closed form in quad precision, where 1 - p and p - q are exact for any long double p
    for (long long it = 0; it < stress::scaled(30000); it++){
        int n1 = stress::rand_int(1, 1000000000), n2 = stress::rand_int(1, 1000000000);
        long double p = stress::rand_int(0, 1) ? 0.5L + stress::rand_int(-(1 << 20), 1 << 20) * powl(10, -stress::rand_int(7, 25)) : stress::rand_int(1, 999999) / 1e6L;
        long double got = gamblers_ruin(n1, n2, p);
        __float128 expected = quad_reference(n1, n2, p);
        assert(0 <= got && got <= 1 && fabsq(got - expected) <= (__float128)1e-12L * expected + (__float128)1e-4900L);
    }

    assert(gamblers_ruin(5, 7, 0.0L) == 1 && gamblers_ruin(5, 7, 1.0L) == 0);
    assert(gamblers_ruin(0, 7, 0.3L) == 1 && gamblers_ruin(5, 0, 0.7L) == 0);
    assert(fabsl(gamblers_ruin(3, 9, 0.5L) - 0.75L) < 1e-18L);
    return 0;
}
