#include "common.h"

#define main library_main
#include "../code_library/fraction.cpp"
#undef main

/// Arithmetic and ordering against exact cross multiplication identities in __int128
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

        Fraction s = a + b, d = a - b, p = a * b;
        __int128 den = (__int128)ad * bd;
        assert((__int128)s.num * den == ((__int128)an * bd + (__int128)bn * ad) * s.den);
        assert((__int128)d.num * den == ((__int128)an * bd - (__int128)bn * ad) * d.den);
        assert((__int128)p.num * den == (__int128)an * bn * p.den);
        for (Fraction f : {s, d, p}) assert(f.den > 0 && __gcd(llabs(f.num), f.den) == 1);
        if (bn != 0){
            Fraction q = a / b;
            assert(q * b == a);
        }
        assert(s - b == a && d + b == a);
        assert((a < b) == ((__int128)an * bd * (ad > 0 ? 1 : -1) * (bd > 0 ? 1 : -1) < (__int128)bn * ad * (ad > 0 ? 1 : -1) * (bd > 0 ? 1 : -1)));
        assert((a == b) == ((__int128)an * bd == (__int128)bn * ad));
        assert((a < b) + (a == b) + (a > b) == 1);
    }

    /// Shared denominators near 1e18: the unreduced denominator is ~1e36, only __int128 holds it
    for (long long it = 0; it < stress::scaled(20000); it++){
        long long den = stress::rand_int(900000000000000000LL, 1000000000000000000LL);
        long long x = stress::rand_int(-den, den), y = stress::rand_int(-den, den);
        Fraction a(x, den), b(y, den), s = a + b, d = a - b;
        assert((__int128)s.num * den == (__int128)(x + y) * s.den && s.den > 0 && __gcd(llabs(s.num), s.den) == 1);
        assert((__int128)d.num * den == (__int128)(x - y) * d.den);
        assert((a < b) == (x < y) && (a == b) == (x == y));
        if (y != 0) assert(a / b == Fraction(x, y));
        long long k = stress::rand_int(-5, 5);
        assert(k + a == a + k && k * a == a * k && (k - a) + a == Fraction(k));
    }
    return 0;
}
