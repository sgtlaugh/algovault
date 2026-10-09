#include "common.h"

#define main library_main
#include "../code_library/integer_root.cpp"
#undef main

/// Floor k-th root by binary search with __int128 powers, no floating point
unsigned long long brute(unsigned long long n, int k){
    if (k == 1) return n;
    unsigned long long lo = 0, hi = n < 2 ? n : min<unsigned long long>(n, 4294967296ULL);
    auto fits = [&](unsigned long long r){
        __int128 p = 1;
        for (int i = 0; i < k; i++){
            p *= r;
            if (p > (__int128)n) return false;
        }
        return true;
    };
    while (lo < hi){
        unsigned long long mid = lo + (hi - lo + 1) / 2;
        if (fits(mid)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

int main(){
    for (unsigned long long n = 0; n <= 100000; n++){
        assert(isqrt(n) == brute(n, 2) && icbrt(n) == brute(n, 3));
    }

    /// Exact powers and their neighbours, where floating estimates go wrong
    for (int k = 1; k <= 64; k++){
        for (unsigned long long r = 1; r <= 3000000; r = r < 100 ? r + 1 : r * 1.37){
            unsigned long long p = 1;
            bool overflow = false;
            for (int i = 0; i < k && !overflow; i++){
                if (p > ULLONG_MAX / r) overflow = true;
                else p *= r;
            }
            if (overflow) break;
            for (unsigned long long n : {p - 1, p, p + 1}) assert(iroot(n, k) == brute(n, k));
        }
    }

    for (long long it = 0; it < stress::scaled(300000); it++){
        unsigned long long n = (unsigned long long)stress::rng()();
        if (it % 3 == 0) n >>= stress::rand_int(0, 63);
        if (it % 5 == 0) n = ULLONG_MAX - stress::rand_int(0, 1000);
        int k = stress::rand_int(1, 64);
        assert(iroot(n, k) == brute(n, k));
        assert(isqrt(n) == brute(n, 2) && icbrt(n) == brute(n, 3));
    }
    return 0;
}
