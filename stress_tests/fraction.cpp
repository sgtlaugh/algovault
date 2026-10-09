#include "common.h"

#define main library_main
#include "../code_library/fraction.cpp"
#undef main

/// Arithmetic and ordering against long double approximations and cross multiplication identities
int main(){
    for (long long it = 0; it < stress::scaled(200000); it++){
        long long range = it % 2 ? 20 : 1000000;
        long long an = stress::rand_int(-range, range), ad = stress::rand_int(1, range) * (stress::rand_int(0, 1) ? 1 : -1);
        long long bn = stress::rand_int(-range, range), bd = stress::rand_int(1, range) * (stress::rand_int(0, 1) ? 1 : -1);
        Fraction a(an, ad), b(bn, bd);

        for (Fraction f : {a, b}){
            assert(f.den > 0 && __gcd(llabs(f.num), f.den) == 1);
        }
        assert((__int128)a.num * ad == (__int128)an * a.den);

        long double x = (long double)an / ad, y = (long double)bn / bd;
        Fraction s = a + b, d = a - b, p = a * b;
        assert(fabsl((long double)s.num / s.den - (x + y)) < 1e-9L);
        assert(fabsl((long double)d.num / d.den - (x - y)) < 1e-9L);
        assert(fabsl((long double)p.num / p.den - x * y) < 1e-9L);
        if (bn != 0){
            Fraction q = a / b;
            assert(q * b == a);
        }
        assert(s - b == a && d + b == a);
        assert((a < b) == ((__int128)an * bd * (ad > 0 ? 1 : -1) * (bd > 0 ? 1 : -1) < (__int128)bn * ad * (ad > 0 ? 1 : -1) * (bd > 0 ? 1 : -1)));
        assert((a == b) == ((__int128)an * bd == (__int128)bn * ad));
        assert((a < b) + (a == b) + (a > b) == 1);
    }
    return 0;
}
