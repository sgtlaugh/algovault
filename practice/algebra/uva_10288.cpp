/***
 *
 * Fraction
 * Exact rational arithmetic, always kept reduced with a positive denominator
 *
 * Complexity: O(log) per operation for the gcd
 *
 * Fraction f(num, den): den != 0, Fraction(x) = x / 1, -1/2 and 1/-2 become the same value
 * f + g, f - g, f * g, f / g (g != 0), -f, +=, -=, *=, /=, comparisons ==, !=, <, <=, >, >=, f.num and f.den
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

        /// 64-bit gcd is about 2x faster, -LLONG_MAX keeps abs(LLONG_MIN) out of std::gcd
        __int128 g = (n >= -LLONG_MAX && n <= LLONG_MAX && d <= LLONG_MAX) ? __gcd(abs((long long)n), (long long)d) : gcd128(n, d);
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

    /// Non-members so an integer works on either side, 1 + half as well as half + 1
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
};

int main(){
    int n;
    while (scanf("%d", &n) == 1){
        Fraction res;
        for (int i = 1; i <= n; i++) res += Fraction(n, i);

        if (res.den == 1){
            printf("%lld\n", res.num);
            continue;
        }

        string whole = to_string(res.num / res.den), rem = to_string(res.num % res.den), den = to_string(res.den);
        string pad(whole.size() + 1, ' ');
        printf("%s%s\n", pad.c_str(), rem.c_str());
        printf("%s %s\n", whole.c_str(), string(max(rem.size(), den.size()), '-').c_str());
        printf("%s%s\n", pad.c_str(), den.c_str());
    }
    return 0;
}
