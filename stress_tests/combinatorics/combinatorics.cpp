#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/combinatorics.cpp"
#undef main

bool is_prime(long long p){
    if (p < 2) return false;
    for (long long d = 2; d * d <= p; d++){
        if (p % d == 0) return false;
    }
    return true;
}

long long random_prime(long long lo, long long hi){
    for (;;){
        long long p = stress::rand_int(lo, hi);
        if (is_prime(p)) return p;
    }
}

/// nCr / nPr / inverses against Pascal's triangle and direct products, for primes just above n up to 2^31 - 1
void check_tables(int n, long long mod){
    Combinatorics comb(n, mod);
    vector<long long> row(1, 1);
    for (int i = 0; i <= n; i++){
        for (int k = -1; k <= i + 1; k++){
            long long pascal = (k < 0 || k > i) ? 0 : row[k];
            assert(comb.nCr(i, k) == pascal);

            long long perm = (k < 0 || k > i) ? 0 : 1;
            for (int j = 0; k >= 0 && j < k && k <= i; j++) perm = perm * (i - j) % mod;
            assert(comb.nPr(i, k) == perm);
        }
        if (i >= 1) assert((__int128)comb.inverse(i) * i % mod == 1);
        assert((__int128)comb.factorial(i) * comb.inv_factorial(i) % mod == 1);

        vector<long long> next(i + 2, 1);
        for (int k = 1; k <= i; k++) next[k] = (row[k - 1] + row[k]) % mod;
        row = next;
    }
}

void check_inverse(long long a, long long m){
    long long inv = mod_inverse(a, m);
    long long r = (long long)(((__int128)a % m + m) % m);

    if (__gcd(r, m) != 1){
        assert(inv == -1);
        return;
    }
    assert(0 <= inv && inv < m);
    assert((__int128)r * inv % m == (m == 1 ? 0 : 1));
}

void check_diophantine(long long a, long long b, long long c){
    long long x, y;
    long long g = __gcd(llabs(a), llabs(b));
    bool solvable = g == 0 ? c == 0 : c % g == 0;
    assert(diophantine(a, b, c, x, y) == solvable);
    if (solvable) assert((__int128)a * x + (__int128)b * y == c);
}

/// reference for wide boxes: scans x and solves y exactly, so only the x range must be small
long long scan_count(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    long long total = 0;
    for (long long x = x1; x <= x2; x++){
        long long rest = c - a * x;
        if (b == 0){
            if (rest == 0) total += max(0LL, y2 - y1 + 1);
            continue;
        }
        if (rest % b == 0 && y1 <= rest / b && rest / b <= y2) total++;
    }
    return total;
}

/// count_diophantine against a double loop over the whole box
void check_count_box(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    long long brute = 0;
    for (long long x = x1; x <= x2; x++){
        for (long long y = y1; y <= y2; y++) brute += a * x + b * y == c;
    }
    assert(count_diophantine(a, b, c, x1, x2, y1, y2) == brute);
}

int main(){
    const long long INT31 = (1LL << 31) - 1, E9 = 1000000000LL;

    for (int n = 0; n <= 40; n++){
        check_tables(n, random_prime(n + 1, 4 * n + 10));
        check_tables(n, INT31);
    }
    for (long long it = 0; it < stress::scaled(20); it++){
        check_tables(stress::rand_int(50, 400), random_prime(401, INT31));
    }

    for (long long m = 1; m <= 60; m++){
        for (long long a = -70; a <= 70; a++) check_inverse(a, m);
    }
    for (long long it = 0; it < stress::scaled(200000); it++){
        long long m = stress::rand_int(0, 3) ? stress::rand_int(1, LLONG_MAX) : stress::rand_int(LLONG_MAX - 1000, LLONG_MAX);
        check_inverse(stress::rand_int(LLONG_MIN, LLONG_MAX), m);
    }

    for (long long a = -20; a <= 20; a++){
        for (long long b = -20; b <= 20; b++){
            for (long long c = -25; c <= 25; c++) check_diophantine(a, b, c);
        }
    }
    for (long long it = 0; it < stress::scaled(200000); it++){
        long long a = stress::rand_int(-E9, E9), b = stress::rand_int(-E9, E9), c = stress::rand_int(-E9, E9);
        if (it % 4 == 0) a = stress::rand_int(0, 1) ? E9 : -E9, b = stress::rand_int(0, 1) ? E9 - 1 : 1 - E9;
        check_diophantine(a, b, c);
    }

    for (long long a = -6; a <= 6; a++){
        for (long long b = -6; b <= 6; b++){
            for (long long c = -12; c <= 12; c++) check_count_box(a, b, c, -4, 5, -5, 3);
        }
    }
    for (long long a = -2; a <= 2; a++){
        for (long long b = -2; b <= 2; b++){
            for (long long c = -3; c <= 3; c++){
                for (long long x1 = -3; x1 <= 3; x1++){
                    for (long long x2 = -3; x2 <= 3; x2++){
                        for (long long y1 = -3; y1 <= 3; y1++){
                            for (long long y2 = -3; y2 <= 3; y2++) check_count_box(a, b, c, x1, x2, y1, y2);
                        }
                    }
                }
            }
        }
    }
    for (long long it = 0; it < stress::scaled(30000); it++){
        long long a = stress::rand_int(-30, 30), b = stress::rand_int(-30, 30), c = stress::rand_int(-200, 200);
        long long x1 = stress::rand_int(-40, 40), x2 = stress::rand_int(-40, 40), y1 = stress::rand_int(-40, 40), y2 = stress::rand_int(-40, 40);
        check_count_box(a, b, c, x1, x2, y1, y2);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long a = stress::rand_int(-E9, E9), b = stress::rand_int(-E9, E9), c = stress::rand_int(-E9, E9);
        if (it % 3 == 0) a = stress::rand_int(-50, 50), b = stress::rand_int(-50, 50);
        if (it % 7 == 0) a = 0;
        if (it % 11 == 0) b = 0;
        if (it % 13 == 0) a = stress::rand_int(0, 1) ? E9 : -E9, b = stress::rand_int(0, 1) ? E9 - 1 : 1 - E9;

        long long x1 = stress::rand_int(-E9, E9 - 2000), x2 = x1 + stress::rand_int(-1, 2000);
        if (it % 5 == 0) x1 = -E9, x2 = -E9 + 2000;
        if (it % 5 == 1) x1 = E9 - 2000, x2 = E9;
        long long y1 = stress::rand_int(0, 1) ? -E9 : stress::rand_int(-E9, E9), y2 = stress::rand_int(0, 1) ? E9 : stress::rand_int(y1 - 1, E9);

        assert(count_diophantine(a, b, c, x1, x2, y1, y2) == scan_count(a, b, c, x1, x2, y1, y2));
        assert(count_diophantine(b, a, c, y1, y2, x1, x2) == scan_count(a, b, c, x1, x2, y1, y2));
    }
    assert(count_diophantine(0, 0, 0, -E9, E9, -E9, E9) == (2 * E9 + 1) * (2 * E9 + 1));

    return 0;
}
