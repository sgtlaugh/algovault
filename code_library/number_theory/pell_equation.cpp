/***
 *
 * Pell Equation
 * Integer solutions of x^2 - d y^2 = 1 and x^2 - d y^2 = -1, exact or modulo m
 *
 * Complexity: O(r + log k) operations on T, for r = O(sqrt(d) log d) the period of the continued fraction of sqrt(d)
 *
 * negative_pell<T>(d): smallest solution (x, y) of x^2 - d y^2 = -1 with y > 0, nullopt if there is none
 * pell<T>(d): smallest solution (x1, y1) of x^2 - d y^2 = 1 with y1 > 0
 * pell<T>(d, k): (x_k, y_k) with x_k + y_k sqrt(d) = (x1 + y1 sqrt(d))^k, for k >= 0
 * pell(d, k, m): (x_k mod m, y_k mod m) as long long, for k >= 0 and 1 <= m < 2^63
 *
 * Notes:
 *   - T is long long by default, any integer type with +, *, % and construction from long long works, like Bignum
 *   - Exact while the answer fits T: x1 fits long long for every d <= 276, __int128 for every d <= 1020
 *   - The (x_k, y_k) for k >= 1 are all the solutions of x^2 - d y^2 = 1 with x, y > 0, in increasing order
 *   - The solutions of x^2 - d y^2 = -1 are the odd powers of negative_pell(d), there are some iff r is odd
 *   - d a perfect square: pell gives (1, 0) for every k, negative_pell gives nullopt, except (0, 1) for d = 1
 *   - 1 <= d < 2^63, but r bounds the speed and x1 has O(r log d) digits
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Exact when m = 0
template<typename T>
T pell_reduce(const T& x, long long m){
    return m ? T(x % m) : x;
}

/// Last convergent p / q of the first period of sqrt(d), p^2 - d q^2 is -1 if the period is odd, else 1
template<typename T>
tuple<T, T, bool> pell_convergent(long long d, long long m){
    assert(d >= 1);
    long long a0 = sqrt((double)d);
    while (a0 * a0 > d) a0--;
    while ((__int128)(a0 + 1) * (a0 + 1) <= d) a0++;
    if (a0 * a0 == d) return {pell_reduce(T(1), m), T(0), false};

    T p0 = pell_reduce(T(1), m), q0 = T(0), p = pell_reduce(T(a0), m), q = pell_reduce(T(1), m);
    long long num = 0, den = 1, a = a0;
    bool odd = false;
    while (true){
        num = den * a - num;
        den = (d - num * num) / den;
        a = (a0 + num) / den;
        odd = !odd;
        if (a == 2 * a0) break;

        T p1 = pell_reduce(a * p + p0, m), q1 = pell_reduce(a * q + q0, m);
        p0 = p, q0 = q, p = p1, q = q1;
    }

    return {p, q, odd};
}

template<typename T = long long>
optional<pair<T, T>> negative_pell(long long d){
    if (d == 1) return pair<T, T>(T(0), T(1));

    auto [p, q, odd] = pell_convergent<T>(d, 0);
    if (!odd) return nullopt;
    return pair<T, T>(p, q);
}

/// With m < 2^63 every value below is under 2^127, so unsigned __int128 never overflows
template<typename T>
pair<T, T> pell_multiply(const pair<T, T>& u, const pair<T, T>& v, long long d, long long m){
    T dyy = pell_reduce(pell_reduce(d * u.second, m) * v.second, m);
    return {pell_reduce(u.first * v.first + dyy, m), pell_reduce(u.first * v.second + u.second * v.first, m)};
}

template<typename T>
pair<T, T> pell_power(long long d, long long k, long long m){
    assert(k >= 0);
    auto [p, q, odd] = pell_convergent<T>(d, m);
    pair<T, T> base(p, q), res(pell_reduce(T(1), m), T(0));
    if (odd) base = pell_multiply(base, base, d, m);

    /// The last squaring is skipped, it would overflow an exact T even when the answer fits
    while (k){
        if (k & 1) res = pell_multiply(res, base, d, m);
        k >>= 1;
        if (k) base = pell_multiply(base, base, d, m);
    }

    return res;
}

template<typename T = long long>
pair<T, T> pell(long long d, long long k = 1){
    return pell_power<T>(d, k, 0);
}

pair<long long, long long> pell(long long d, long long k, long long m){
    assert(m >= 1);
    auto [x, y] = pell_power<unsigned __int128>(d, k, m);  /// unsigned measured 9% faster
    return {(long long)x, (long long)y};
}

int main(){
    /// 3^2 - 2 * 2^2 = 1, (3 + 2 sqrt(2))^2 = 17 + 12 sqrt(2) and (3 + 2 sqrt(2))^3 = 99 + 70 sqrt(2)
    assert((pell(2) == make_pair(3LL, 2LL)));
    assert((pell(2, 2) == make_pair(17LL, 12LL)));
    assert((pell(2, 3, 50) == make_pair(49LL, 20LL)));

    assert((pell(61) == make_pair(1766319049LL, 226153980LL)));
    assert((pell<__int128>(61) == make_pair((__int128)1766319049, (__int128)226153980)));

    /// sqrt(13) has period 5, odd, so 18^2 - 13 * 5^2 = -1 and pell(13) = (18 + 5 sqrt(13))^2 = 649 + 180 sqrt(13)
    assert((negative_pell(13) == make_pair(18LL, 5LL)));
    assert((pell(13) == make_pair(649LL, 180LL)));
    assert(!negative_pell(3));  /// sqrt(3) has period 2, x^2 = 3 y^2 - 1 is 2 mod 3, never a square

    return 0;
}
