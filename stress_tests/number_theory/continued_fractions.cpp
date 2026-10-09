#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/continued_fractions.cpp"
#undef main

using Frac = pair<long long, long long>;

const long long E18 = 1000000000000000000LL;

/// Rebuilds the value from the back of the expansion, independent of the convergent recurrence
Frac evaluate(const vector<long long>& a){
    __int128 p = a.back(), q = 1;
    for (int i = (int)a.size() - 2; i >= 0; i--){
        swap(p, q);
        p += (__int128)a[i] * q;
    }
    return {(long long)p, (long long)q};
}

Frac reduce(long long p, long long q){
    long long g = __gcd(p < 0 ? -p : p, q);
    return {p / g, q / g};
}

long long floor_div(__int128 a, __int128 b){
    __int128 r = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) r--;
    return (long long)r;
}

/// Scans denominators and numerators upwards keeping only strict improvements, so ties go to the smaller denominator, then the smaller value
Frac brute_best(long long p, long long q, long long n){
    Frac best = {0, 0};
    __int128 best_d = 0;

    for (long long b = 1; b <= n; b++){
        long long lo = floor_div((__int128)p * b, q);
        for (long long a : {lo, lo + 1}){
            __int128 d = (__int128)p * b - (__int128)a * q;
            if (d < 0) d = -d;
            if (best.second == 0){
                best = {a, b}, best_d = d;
                continue;
            }

            if (d * best.second < best_d * b) best = {a, b}, best_d = d;
        }
    }

    return best;
}

/// Smallest fraction with p, q <= n satisfying f, over every pair
template<typename F>
Frac brute_search(F f, long long n){
    Frac best = {1, 0};
    for (long long q = 1; q <= n; q++){
        for (long long p = 0; p <= n; p++){
            if (f(p, q) && (best.second == 0 || p * best.second < best.first * q)) best = {p, q};
        }
    }
    return best.second ? reduce(best.first, best.second) : best;
}

/// Smallest fraction >= tp / tq with denominator <= n, for 0 <= tp / tq < 1 with reduced denominator > n
/// One of the last convergent within n and the largest semiconvergent lies above tp / tq, that one is the answer
Frac right_neighbour(long long tp, long long tq, long long n){
    auto c = convergents(continued_fraction(tp, tq));
    size_t k = 0;
    while (k + 1 < c.size() && c[k + 1].second <= n) k++;

    Frac prev = k ? c[k - 1] : Frac{1, 0};
    long long t = (n - prev.second) / c[k].second;
    Frac semi = {prev.first + t * c[k].first, prev.second + t * c[k].second};
    bool conv_above = (__int128)c[k].first * tq > (__int128)tp * c[k].second;
    return conv_above ? c[k] : semi;
}

/// The header's 4 log2(N) + 10, with log2 rounded down
long long call_bound(long long n){
    return 4 * (63 - __builtin_clzll(n)) + 10;
}

void check_expansion(long long p, long long q){
    auto a = continued_fraction(p, q);
    auto c = convergents(a);
    assert(evaluate(a) == reduce(p, q) && c.back() == reduce(p, q));
    assert(a[0] == floor_div(p, q));

    for (size_t i = 1; i < a.size(); i++) assert(a[i] >= 1);
    if (a.size() > 1) assert(a.back() >= 2);

    for (size_t i = 0; i < c.size(); i++){
        assert(c[i] == evaluate(vector<long long>(a.begin(), a.begin() + i + 1)));
        if (i) assert((__int128)c[i].first * c[i - 1].second - (__int128)c[i - 1].first * c[i].second == (i % 2 ? 1 : -1));
    }
}

int main(){
    for (long long q = 1; q <= 60; q++){
        for (long long p = -150; p <= 150; p++) check_expansion(p, q);
    }

    for (long long it = 0; it < stress::scaled(20000); it++){
        long long q = stress::rand_int(1, it % 2 ? E18 : 1000000);
        long long p = stress::rand_int(-E18, E18);
        check_expansion(p, q);
    }

    for (long long it = 0; it < stress::scaled(2000); it++){
        long long q = stress::rand_int(1, it % 2 ? LLONG_MAX : 1000);
        check_expansion(LLONG_MIN + 1 + stress::rand_int(0, 5), q);
        check_expansion(LLONG_MAX - stress::rand_int(0, 5), q);
    }

    for (long long q = 1; q <= 25; q++){
        for (long long p = -60; p <= 60; p++){
            for (long long n = 1; n <= 30; n++) assert(best_approximation(p, q, n) == brute_best(p, q, n));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long q = stress::rand_int(1, it % 3 ? E18 : 100000);
        long long p = stress::rand_int(-E18, E18);
        long long n = stress::rand_int(1, it % 2 ? 300 : 3000);
        assert(best_approximation(p, q, n) == brute_best(p, q, n));
    }

    /// Exhaustive search over thresholds tp / tq with >= and >, and irrational thresholds sqrt(c)
    for (long long n = 1; n <= 14; n++){
        for (long long tq = 1; tq <= 16; tq++){
            for (long long tp = 0; tp <= 18; tp++){
                auto ge = [&](long long p, long long q){ return p * tq >= tp * q; };
                auto gt = [&](long long p, long long q){ return p * tq > tp * q; };
                assert(fraction_search(ge, n) == brute_search(ge, n));
                assert(fraction_search(gt, n) == brute_search(gt, n));
            }
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long n = stress::rand_int(1, it % 2 ? 40 : 120);
        long long c = stress::rand_int(1, 200);
        long long calls = 0;
        auto f = [&](long long p, long long q){
            calls++;
            assert(p >= 0 && p <= n && q >= 1 && q <= n && __gcd(p, q) == 1);
            return p * p >= c * q * q;
        };
        assert(fraction_search(f, n) == brute_search([&](long long p, long long q){ return p * p >= c * q * q; }, n));
        assert(calls <= call_bound(n));
    }

    /// Partial quotients all 1 (golden ratio), all 2 (1 + sqrt 2) and alternating 1, 2 ((1 + sqrt 3) / 2) cost the most calls per bit of N
    auto phi = [](long long p, long long q){ return (__int128)p * p - (__int128)p * q - (__int128)q * q >= 0; };
    auto silver = [](long long p, long long q){ return (__int128)p * p - (__int128)2 * p * q - (__int128)q * q >= 0; };
    auto half_sqrt3 = [](long long p, long long q){ return (__int128)(2 * p - q) * (2 * p - q) >= (__int128)3 * q * q; };
    auto count_calls = [](auto f, long long n){
        long long calls = 0;
        fraction_search([&](long long p, long long q){ calls++; return f(p, q); }, n);
        return calls;
    };

    for (long long n = 1; n <= E18; n += n / 64 + stress::rand_int(1, 3)){
        assert(count_calls(phi, n) <= call_bound(n));
        assert(count_calls(silver, n) <= call_bound(n));
        assert(count_calls(half_sqrt3, n) <= call_bound(n));
    }

    /// Large bounds: the Stern-Brocot search and the convergent construction must agree
    for (long long it = 0; it < stress::scaled(20000); it++){
        long long n = it % 4 ? stress::rand_int(1, E18) : E18 - stress::rand_int(0, 10);
        long long tq = stress::rand_int(2, E18), tp = stress::rand_int(0, tq - 1);
        if (it % 5 == 0) tq = stress::rand_int(2, 1000), tp = stress::rand_int(0, tq - 1);
        long long calls = 0;
        auto ge = [&](long long p, long long q){
            calls++;
            assert(p >= 0 && p <= n && q >= 1 && q <= n);
            return (__int128)p * tq >= (__int128)tp * q;
        };

        Frac got = fraction_search(ge, n), t = reduce(tp, tq);
        assert(got == (t.second <= n ? t : right_neighbour(tp, tq, n)));
        assert(calls <= call_bound(n));

        /// The best approximation is the closer of the two neighbours, the left one found by mirroring around 1 / 2
        if (t.second <= n) continue;
        auto le_mirror = [&](long long p, long long q){ return (__int128)p * tq >= (__int128)(tq - tp) * q; };
        Frac r = got, mirrored = fraction_search(le_mirror, n), l = {mirrored.second - mirrored.first, mirrored.second};

        __int128 dl = (__int128)tp * l.second - (__int128)l.first * tq, dr = (__int128)r.first * tq - (__int128)tp * r.second;
        __int128 lhs = dl * r.second, rhs = dr * l.second;
        Frac expected = lhs < rhs || (lhs == rhs && l.second <= r.second) ? l : r;
        assert(best_approximation(tp, tq, n) == expected);

        long long shift = stress::rand_int(-(E18 - tq) / tq, (E18 - tq) / tq);
        assert(best_approximation(tp + shift * tq, tq, n) == Frac(expected.first + shift * expected.second, expected.second));
    }

    return 0;
}
