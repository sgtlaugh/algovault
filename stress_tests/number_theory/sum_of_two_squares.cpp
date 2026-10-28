#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/sum_of_two_squares.cpp"
#undef main

using Pairs = vector<pair<long long, long long>>;

long long isqrt(long long n){
    long long r = sqrtl(n);
    while (r * r > n) r--;
    while ((__int128)(r + 1) * (r + 1) <= n) r++;
    return r;
}

Pairs brute(long long n){
    Pairs res;
    for (long long x = 0; x * x <= n; x++){
        long long y = isqrt(n - x * x);
        if (x * x + y * y == n) res.push_back({x, y});
    }
    return res;
}

/// Pairs must solve n, be sorted without repeats, and number r2 / 4 plus one more (0, x) when n is a square
void check_valid(long long n, const Pairs& got, long long r2){
    for (auto [x, y]: got){
        assert(x >= 0 && y >= 0);
        assert((__int128)x * x + (__int128)y * y == n);
    }
    for (size_t i = 0; i + 1 < got.size(); i++) assert(got[i] < got[i + 1]);

    long long r = isqrt(n);
    assert((long long)got.size() == r2 / 4 + (r * r == n));
}

long long random_prime(long long lo, long long hi, int residue){
    while (true){
        long long p = stress::rand_int(lo, hi);
        if (p % 4 == residue && prm::is_prime(p)) return p;
    }
}

int main(){
    TwoSquares squares;

    /// Every n up to N against all (x, y) bucketed by x^2 + y^2
    const int N = 200000;
    vector<Pairs> all(N + 1);
    for (long long x = 0; x * x <= N; x++){
        for (long long y = 0; x * x + y * y <= N; y++) all[x * x + y * y].push_back({x, y});
    }
    for (int n = 0; n <= N; n++){
        assert(squares.pairs(n) == all[n]);
        long long r2 = 0;
        for (auto [x, y]: all[n]) r2 += (x == 0 ? 1 : 2) * (y == 0 ? 1 : 2);
        assert(squares.r2(n) == r2);
    }

    /// Random n up to 10^12, where large prime factors and Pollard rho splits appear
    for (int it = 0; it < stress::scaled(20); it++){
        long long n = stress::rand_int(1, 1000000000000LL);
        auto got = squares.pairs(n);
        assert(got == brute(n));
        check_valid(n, got, squares.r2(n));
    }

    /// Known factorizations up to 2^63: r2 = 4 prod (e + 1) over primes 1 mod 4, 0 if a prime 3 mod 4 has odd e
    for (int it = 0; it < stress::scaled(3000); it++){
        __int128 n = 1;
        long long r2 = 4;
        map<long long, int> exponent;
        while (true){
            int shape = stress::rand_int(0, 5);
            long long p;
            if (shape == 0) p = 2;
            else if (shape == 1) p = random_prime(3, 100, 3);
            else if (shape == 2) p = random_prime(3, 1000000000, 3);
            else if (shape == 3) p = random_prime(5, 100, 1);
            else if (shape == 4) p = random_prime(5, 1000000000, 1);
            else p = random_prime(5, 4000000000000000000LL, 1);

            int e = stress::rand_int(1, 3);
            __int128 m = n;
            for (int k = 0; k < e && m <= LLONG_MAX; k++) m *= p;
            if (m > LLONG_MAX) break;

            n = m;
            exponent[p] += e;
        }

        for (auto [p, e]: exponent){
            if (p % 4 == 1) r2 *= e + 1;
            if (p % 4 == 3 && e % 2 == 1) r2 = 0;
        }
        assert(squares.r2(n) == r2);
        check_valid(n, squares.pairs(n), r2);
    }

    /// Edge cases: 2, 3, large primes 1 and 3 mod 4, 2^62, the 11 smallest primes 1 mod 4, the most pairs below 2^63
    assert(squares.pairs(3).empty());
    assert(squares.pairs(2) == Pairs({{1, 1}}));
    long long p3 = random_prime(4000000000000000000LL, LLONG_MAX, 3), p1 = random_prime(4000000000000000000LL, LLONG_MAX, 1);
    assert(squares.pairs(p3).empty());
    check_valid(p1, squares.pairs(p1), 8);
    check_valid(1LL << 62, squares.pairs(1LL << 62), 4);

    long long most = 1;
    for (long long p: {5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97}) most *= p;
    for (long long k: {1LL, 5LL, 9LL}) check_valid(most * k, squares.pairs(most * k), squares.r2(most * k));
    assert(squares.r2(most) == 4 << 11);

    long long most_63 = 6437978455413290825LL;
    check_valid(most_63, squares.pairs(most_63), squares.r2(most_63));
    assert(squares.pairs(most_63).size() == 4608);

    return 0;
}
