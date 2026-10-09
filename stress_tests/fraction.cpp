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

    /// Independent numerators and denominators up to ~1e18 with a large shared factor, so sums still fit after reducing
    for (long long it = 0; it < stress::scaled(20000); it++){
        long long g = stress::rand_int(100000000000000000LL, 300000000000000000LL);
        long long ad = g * stress::rand_int(1, 3), bd = g * stress::rand_int(1, 3);
        long long an = stress::rand_int(-ad, ad), bn = stress::rand_int(-bd, bd);
        Fraction a(an, ad), b(bn, bd), s = a + b, d = a - b;
        long long lcm = ad / __gcd(ad, bd) * bd;
        __int128 x = (__int128)an * (lcm / ad), y = (__int128)bn * (lcm / bd);
        assert((__int128)s.num * lcm == (x + y) * s.den && (__int128)d.num * lcm == (x - y) * d.den);
        for (Fraction f : {a, b, s, d}) assert(f.den > 0 && __gcd(llabs(f.num), f.den) == 1);
        assert(s - b == a && d + b == a);
    }

    Fraction low(LLONG_MIN, 2);
    assert(low.num == LLONG_MIN / 2 && low.den == 1 && Fraction(LLONG_MIN, 1).num == LLONG_MIN);
    /// The unreduced numerator 3 * 2^62 overflows long long, the 64-bit gcd path must not see it
    assert(Fraction(1LL << 62, 3) * Fraction(3) == Fraction(1LL << 62) && Fraction(-(1LL << 62), 3) * 3 == Fraction(-(1LL << 62)));
    return 0;
}
