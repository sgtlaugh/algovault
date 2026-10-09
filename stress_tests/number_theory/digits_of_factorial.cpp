// LINK: -lgmp -lquadmath
#include "../common.h"
#include <gmp.h>
#include <quadmath.h>

#define main library_main
#include "../../code_library/number_theory/digits_of_factorial.cpp"
#undef main

long long exact_digits(const mpz_t f, long long b){
    mpz_t p;
    mpz_init_set_ui(p, 1);
    long long d = 0;
    while (mpz_cmp(p, f) <= 0) mpz_mul_ui(p, p, b), d++;
    mpz_clear(p);
    return d;
}

/// Beyond n = 20 the library rounds a long double logarithm, so a log within a few ulp of an integer is out of its reach
bool decidable(long long n, long long b){
    __float128 l = lgammaq((__float128)n + 1) / logq(b);
    return n <= 20 || fabsq(l - roundq(l)) > (__float128)1e-17 * l + (__float128)1e-14;
}

int main(){
    long long skipped = 0;
    mpz_t f, r;
    mpz_inits(f, r, NULL);

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long n = stress::rand_int(0, it % 3 ? 25 : 3000);
        mpz_fac_ui(f, n);

        vector<long long> bases = {stress::rand_int(2, 100), stress::rand_int(2, LLONG_MAX)};
        for (int k = 1; k <= 4; k++){
            mpz_root(r, f, k == 1 ? 1 : stress::rand_int(2, 60));  /// b near a root of n! puts the answer on a knife edge
            if (mpz_cmp_ui(r, 1) > 0 && mpz_sizeinbase(r, 2) < 63){
                long long b = mpz_get_si(r);
                bases.push_back(b);
                if (b > 2) bases.push_back(b - 1);
                if (b < LLONG_MAX) bases.push_back(b + 1);
            }
        }

        for (long long b : bases){
            if (!decidable(n, b)){
                skipped++;
                continue;
            }
            assert(digits_of_factorial(n, b) == exact_digits(f, b));
        }
    }

    auto check_large = [&](long long n, long long b){
        if (!decidable(n, b)) skipped++;
        else assert(digits_of_factorial(n, b) == (long long)floorq(lgammaq((__float128)n + 1) / logq(b)) + 1);
    };

    for (long long it = 0; it < stress::scaled(10000); it++){
        check_large(stress::rand_int(21, 2000000000), stress::rand_int(2, it % 2 ? 1000 : LLONG_MAX));
    }

    /// b just off the k-th root of n! puts log_b(n!) within ~1e-6 of k, past what a double logarithm resolves
    for (long long it = 0; it < stress::scaled(10000); it++){
        long long n = stress::rand_int(1000000, 2000000000);
        __float128 lf = lgammaq((__float128)n + 1), k = floorq(lf / logq(stress::rand_int(1000000000000LL, 4000000000000000000LL)));
        check_large(n, (long long)expq(lf / k) + stress::rand_int(-30000, 30000));
    }

    mpz_clears(f, r, NULL);
    fprintf(stderr, "skipped %lld undecidable cases\n", skipped);
    return 0;
}
