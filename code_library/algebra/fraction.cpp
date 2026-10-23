/***
 *
 * Fraction
 * Exact rational arithmetic in long long, always kept reduced with a positive denominator
 *
 * Complexity: O(log v) per arithmetic operation for the gcd, v the largest numerator or denominator, O(1) per comparison
 *
 * Fraction f(x): the integer x, 0 by default
 * Fraction f(num, den): num / den in lowest terms, den != 0
 * f.num, f.den: the reduced numerator and denominator, den > 0, so 1/-2 and -2/4 both hold -1 and 2
 *
 * f + g, f - g, f * g: the sum, difference and product
 * f / g: the quotient, g != 0
 * -f and f += g, -=, *=, /=: negation and the compound assignments
 * f == g, !=, <, <=, >, >=: exact comparison
 *
 * Every binary operator and comparison also takes a long long on either side
 * Intermediate products use __int128, every reduced result must fit in long long (asserted)
 *
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Fraction{
    long long num, den;

    Fraction(long long x = 0) : num(x), den(1) {}

    Fraction(long long n, long long d){
        *this = make(n, d);
    }

    friend Fraction operator+(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.den + (__int128)b.num * a.den, (__int128)a.den * b.den);
    }

    friend Fraction operator-(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.den - (__int128)b.num * a.den, (__int128)a.den * b.den);
    }

    friend Fraction operator*(const Fraction& a, const Fraction& b){
        return make((__int128)a.num * b.num, (__int128)a.den * b.den);
    }

    friend Fraction operator/(const Fraction& a, const Fraction& b){
        assert(b.num != 0);
        return make((__int128)a.num * b.den, (__int128)a.den * b.num);
    }

    friend Fraction operator-(const Fraction& a){
        return make(-(__int128)a.num, a.den);
    }

    Fraction& operator+=(const Fraction& b){
        return *this = *this + b;
    }

    Fraction& operator-=(const Fraction& b){
        return *this = *this - b;
    }

    Fraction& operator*=(const Fraction& b){
        return *this = *this * b;
    }

    Fraction& operator/=(const Fraction& b){
        return *this = *this / b;
    }

    friend bool operator==(const Fraction& a, const Fraction& b){
        return a.num == b.num && a.den == b.den;
    }

    friend bool operator!=(const Fraction& a, const Fraction& b){
        return !(a == b);
    }

    friend bool operator<(const Fraction& a, const Fraction& b){
        return (__int128)a.num * b.den < (__int128)b.num * a.den;
    }

    friend bool operator>(const Fraction& a, const Fraction& b){
        return b < a;
    }

    friend bool operator<=(const Fraction& a, const Fraction& b){
        return !(b < a);
    }

    friend bool operator>=(const Fraction& a, const Fraction& b){
        return !(a < b);
    }

    static Fraction make(__int128 n, __int128 d){
        assert(d != 0);
        if (d < 0) n = -n, d = -d;

        __int128 g = (n >= -LLONG_MAX && n <= LLONG_MAX && d <= LLONG_MAX) ? gcd((long long)n, (long long)d) : gcd128(n, d);
        n /= g, d /= g;
        assert(n >= LLONG_MIN && n <= LLONG_MAX && d <= LLONG_MAX);

        Fraction f;
        f.num = (long long)n, f.den = (long long)d;
        return f;
    }

    static __int128 gcd128(__int128 a, __int128 b){
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        while (b) a %= b, swap(a, b);
        return a;
    }
};

int main(){
    /// 1/2 and 1/3 under each operation
    Fraction half(1, 2), third(1, 3);
    assert(half + third == Fraction(5, 6));
    assert(half - third == Fraction(1, 6));
    assert(half * third == Fraction(1, 6));
    assert(half / third == Fraction(3, 2));
    assert(third < half);

    Fraction f(6, -8);
    assert(f.num == -3);
    assert(f.den == 4);

    /// An integer works on either side
    assert(1 - half == half);
    assert(2 * third == Fraction(2, 3));
    assert(half < 1);

    Fraction sum;
    for (int k = 1; k <= 4; k++) sum += Fraction(1, k);
    assert(sum == Fraction(25, 12));

    /// The cross products near 10^36 only fit in __int128
    const long long BIG = 1000000000000000000LL;
    assert(Fraction(BIG - 1, BIG) < Fraction(BIG, BIG + 1));

    return 0;
}
