#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/extended_gcd.cpp"
#undef main

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

void check_extended_gcd(long long a, long long b){
    auto [g, x, y] = extended_gcd(a, b);
    assert(g == __gcd(a, b));
    assert((__int128)a * x + (__int128)b * y == g);
}

void check_diophantine(long long a, long long b, long long c){
    long long g = __gcd(llabs(a), llabs(b));
    bool solvable = g == 0 ? c == 0 : c % g == 0;
    auto sol = diophantine(a, b, c);
    assert(sol.has_value() == solvable);
    if (solvable) assert((__int128)a * sol->first + (__int128)b * sol->second == c);
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
    const long long E9 = 1000000000LL;

    for (long long a = 0; a <= 60; a++){
        for (long long b = 0; b <= 60; b++) check_extended_gcd(a, b);
    }
    for (long long it = 0; it < stress::scaled(200000); it++){
        long long a = stress::rand_int(0, LLONG_MAX), b = stress::rand_int(0, it % 2 ? LLONG_MAX : 1000);
        check_extended_gcd(a, b);
        check_extended_gcd(b, a);
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
