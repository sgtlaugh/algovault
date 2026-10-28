#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/pell_equation.cpp"
#undef main

#define main bignum_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

typedef pair<Bignum, Bignum> Big;
typedef pair<long long, long long> Mod;

Bignum big(__int128 v){
    const long long E18 = 1000000000000000000LL;
    Bignum top((long long)(v / E18 / E18)), mid((long long)(v / E18 % E18)), low((long long)(v % E18));
    return (top * E18 + mid) * E18 + low;
}

Big times(const Big& u, const Big& v, long long d){
    return {u.first * v.first + Bignum(d) * u.second * v.second, u.first * v.second + u.second * v.first};
}

Mod times_mod(const Mod& u, const Mod& v, long long d, long long m){
    __int128 dm = d % m;
    __int128 x = ((__int128)u.first * v.first % m + dm * u.second % m * v.second % m) % m;
    __int128 y = ((__int128)u.first * v.second + (__int128)u.second * v.first) % m;
    return {(long long)x, (long long)y};
}

/// The fundamental solution is the first convergent of sqrt(d) with norm 1, found without the period parity shortcut
Big first_unit_convergent(long long d){
    long long a0 = 0;
    while ((a0 + 1) * (a0 + 1) <= d) a0++;
    if (a0 * a0 == d) return {Bignum(1), Bignum(0)};

    Bignum p0 = 1, q0 = 0, p = a0, q = 1;
    long long num = 0, den = 1, a = a0;
    while (p * p - Bignum(d) * q * q != Bignum(1)){
        num = den * a - num;
        den = (d - num * num) / den;
        a = (a0 + num) / den;
        Bignum p1 = Bignum(a) * p + p0, q1 = Bignum(a) * q + q0;
        p0 = p, q0 = q, p = p1, q = q1;
    }

    return {p, q};
}

long long floor_sqrt(long long n){
    long long r = sqrtl(n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

/// Every positive solution with y <= Y by trying each y, against pell(d, k) and the odd powers of negative_pell(d)
void brute_small_d(){
    const long long Y = 30000;
    for (long long d = 1; d <= 300; d++){
        vector<Mod> plus, minus;
        for (long long y = 1; y <= Y; y++){
            long long x = floor_sqrt(d * y * y + 1);
            if (x * x == d * y * y + 1) plus.push_back({x, y});
            x = floor_sqrt(d * y * y - 1);
            if (x * x == d * y * y - 1) minus.push_back({x, y});
        }

        Big fund = pell<Bignum>(d);
        assert(fund == first_unit_convergent(d));
        long long sq = floor_sqrt(d);
        if (sq * sq == d){
            assert(plus.empty() && fund == Big(1, 0));
            assert((pell<Bignum>(d, 5) == Big(1, 0)));
            auto neg = negative_pell<Bignum>(d);
            assert(d == 1 ? neg == Big(0, 1) : !neg);
            continue;
        }

        for (long long k = 1; ; k++){
            Big sol = pell<Bignum>(d, k);
            if (sol.second > Bignum(Y)){
                assert(k == (long long)plus.size() + 1);
                break;
            }
            assert(k <= (long long)plus.size());
            assert(sol == Big(plus[k - 1].first, plus[k - 1].second));
        }

        auto neg = negative_pell<Bignum>(d);
        if (!neg){
            assert(minus.empty());
            continue;
        }

        Big sol = *neg, square = times(*neg, *neg, d);
        assert(times(sol, sol, d) == fund);
        for (size_t j = 0; ; j++){
            if (sol.second > Bignum(Y)){
                assert(j == minus.size());
                break;
            }
            assert(j < minus.size());
            assert(sol == Big(minus[j].first, minus[j].second));
            sol = times(sol, square, d);
        }
    }
}

/// long long holds x1 for d <= 276 and __int128 for d <= 1020, the header's limits
void exact_types(){
    for (long long d = 1; d <= 1020; d++){
        Big ref = first_unit_convergent(d);
        assert(pell<Bignum>(d) == ref);
        assert(ref.first * ref.first - Bignum(d) * ref.second * ref.second == Bignum(1));

        auto [x, y] = pell<__int128>(d);
        assert(big(x) == ref.first && big(y) == ref.second);
        if (d <= 276){
            auto [lx, ly] = pell(d);
            assert(Bignum(lx) == ref.first && Bignum(ly) == ref.second);
        }

        auto neg = negative_pell<__int128>(d);
        auto big_neg = negative_pell<Bignum>(d);
        assert(neg.has_value() == big_neg.has_value());
        if (!neg) continue;
        assert(big(neg->first) == big_neg->first && big(neg->second) == big_neg->second);
        assert(times(*big_neg, *big_neg, d) == ref);
        assert(big_neg->first * big_neg->first - Bignum(d) * big_neg->second * big_neg->second == Bignum(-1));
    }

    assert(pell<Bignum>(277).first > Bignum(LLONG_MAX));
    assert(pell<Bignum>(1021).first > big(numeric_limits<__int128>::max()));
}

/// pell(d, k) against repeated multiplication, pell(d, k, m) against the exact value mod m
void kth_and_mod(){
    vector<long long> mods = {1, 2, 1000000007, LLONG_MAX, LLONG_MAX - 24, 4611686018427387847LL};
    for (long long it = 0; it < stress::scaled(60); it++){
        long long d = it < 10 ? stress::rand_int(1, 30) : stress::rand_int(1, it < 50 ? 3000 : 200000);
        Big fund = pell<Bignum>(d), sol(Bignum(1), Bignum(0));
        if (it % 10 == 0) assert(fund == first_unit_convergent(d));

        vector<long long> ms = mods;
        ms.push_back(stress::rand_int(1, 100));
        ms.push_back(stress::rand_int(1, LLONG_MAX));
        for (long long k = 0; k <= 12; k++){
            assert(pell<Bignum>(d, k) == sol);
            for (long long m : ms){
                Mod got = pell(d, k, m);
                assert(Bignum(got.first) == sol.first % Bignum(m));
                assert(Bignum(got.second) == sol.second % Bignum(m));
            }
            sol = times(sol, fund, d);
        }
    }
}

/// Large k by splitting k = k1 + k2, large d through families with known fundamental solutions
void large_mod(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long d = it % 3 ? stress::rand_int(1, 5000) : stress::rand_int(1, 2000000);
        long long m = it % 4 ? stress::rand_int(1, LLONG_MAX) : stress::rand_int(1, 1000);
        long long k1 = stress::rand_int(0, LLONG_MAX / 2), k2 = stress::rand_int(0, LLONG_MAX / 2);
        if (it % 5 == 0) k1 = stress::rand_int(0, 3);
        assert(pell(d, k1 + k2, m) == times_mod(pell(d, k1, m), pell(d, k2, m), d, m));

        Mod x = pell(d, k1, m);
        __int128 norm = (__int128)x.first * x.first % m - (__int128)(d % m) * x.second % m * x.second % m;
        assert(((norm % m + m) % m) == 1 % m);
    }

    /// a^2 - (a^2 - 1) 1^2 = 1, (2a^2 + 1)^2 - (a^2 + 1)(2a)^2 = 1 and (a^2 + 1)^2 - (a^2 + 2) a^2 = 1
    const long long A = 3037000499LL;
    for (long long it = 0; it < stress::scaled(300); it++){
        long long a = it < 20 ? A - it : stress::rand_int(2, A);
        long long m = it % 2 ? LLONG_MAX - stress::rand_int(0, 100) : stress::rand_int(1, LLONG_MAX);
        long long k = stress::rand_int(0, LLONG_MAX);
        __int128 a2 = (__int128)a * a;
        vector<pair<long long, Mod>> cases = {
            {(long long)(a2 - 1), {a % m, 1 % m}},
            {(long long)(a2 + 1), {(long long)((2 * a2 + 1) % m), (2 * a) % m}},
            {(long long)(a2 + 2), {(long long)((a2 + 1) % m), a % m}},
        };

        for (auto [d, fund] : cases){
            assert(pell(d, 1, m) == fund);
            Mod res(1 % m, 0), base = fund;
            for (long long e = k; e; e >>= 1){
                if (e & 1) res = times_mod(res, base, d, m);
                base = times_mod(base, base, d, m);
            }
            assert(pell(d, k, m) == res);
        }

        long long sq = (long long)a2;
        assert(pell(sq, k, m) == Mod(1 % m, 0));
        assert(pell(sq) == Mod(1, 0));
        assert(!negative_pell(sq));
    }
}

int main(){
    brute_small_d();
    exact_types();
    kth_and_mod();
    large_mod();
    stress::expect_assert("number_theory/pell_equation.cpp", []{ pell(0); });
    stress::expect_assert("number_theory/pell_equation.cpp", []{ pell(2, -1); });
    stress::expect_assert("number_theory/pell_equation.cpp", []{ pell(2, 1, 0); });
    return 0;
}
