/***
 *
 * Fraction
 * Exact rational arithmetic, always kept reduced with a positive denominator
 *
 * Complexity: O(log) per operation for the gcd
 *
 * Fraction f(num, den): den != 0, Fraction(x) = x / 1, -1/2 and 1/-2 become the same value
 * f + g, f - g, f * g, f / g (g != 0), comparisons ==, !=, <, <=, >, >=, f.num and f.den
 * Intermediate products use __int128, the reduced result must fit in long long (asserted)
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Fraction{
    long long num, den;

    static __int128 gcd128(__int128 a, __int128 b){
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        while (b) a %= b, swap(a, b);
        return a;
    }

    static Fraction make(__int128 n, __int128 d){
        assert(d != 0);
        if (d < 0) n = -n, d = -d;
        __int128 g = gcd128(n, d);
        n /= g, d /= g;
        assert(n >= LLONG_MIN && n <= LLONG_MAX && d <= LLONG_MAX);
        Fraction f;
        f.num = (long long)n, f.den = (long long)d;
        return f;
    }

    Fraction(long long x = 0) : num(x), den(1) {}

    Fraction(long long n, long long d){
        *this = make(n, d);
    }

    Fraction operator+(const Fraction& o) const{
        return make((__int128)num * o.den + (__int128)o.num * den, (__int128)den * o.den);
    }

    Fraction operator-(const Fraction& o) const{
        return make((__int128)num * o.den - (__int128)o.num * den, (__int128)den * o.den);
    }

    Fraction operator*(const Fraction& o) const{
        return make((__int128)num * o.num, (__int128)den * o.den);
    }

    Fraction operator/(const Fraction& o) const{
        assert(o.num != 0);
        return make((__int128)num * o.den, (__int128)den * o.num);
    }

    bool operator==(const Fraction& o) const{
        return num == o.num && den == o.den;
    }

    bool operator!=(const Fraction& o) const{
        return !(*this == o);
    }

    bool operator<(const Fraction& o) const{
        return (__int128)num * o.den < (__int128)o.num * den;
    }

    bool operator>(const Fraction& o) const{
        return o < *this;
    }

    bool operator<=(const Fraction& o) const{
        return !(o < *this);
    }

    bool operator>=(const Fraction& o) const{
        return !(*this < o);
    }
};

int main(){
    Fraction half(1, 2), third(1, 3);
    assert(half + third == Fraction(5, 6));
    assert(half - third == Fraction(1, 6));
    assert(half * third == Fraction(1, 6));
    assert(half / third == Fraction(3, 2));
    assert(Fraction(1, -2) == Fraction(-1, 2) && Fraction(-1, 2).den == 2 && Fraction(-1, 2).num == -1);
    assert(Fraction(0, -5) == Fraction(0) && Fraction(0, 7).den == 1);
    assert(Fraction(6, 8).num == 3 && Fraction(6, 8).den == 4);
    assert(third < half && Fraction(-1, 2) < third && !(half < half) && half <= half && half >= third);
    assert(Fraction(4) / Fraction(2) == Fraction(2));

    const long long BIG = 1000000000000000000LL;
    assert(Fraction(BIG, 3) * Fraction(3, BIG) == Fraction(1));
    assert(Fraction(BIG - 1, BIG) < Fraction(BIG, BIG + 1));
    return 0;
}
