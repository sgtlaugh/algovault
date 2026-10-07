/***
 *
 * Shared helpers for the stress tests
 *
 * Each stress test includes this header, then the library file with its main renamed:
 *
 *     #include "common.h"
 *     #define main library_main
 *     #include "../code_library/fenwick_tree.cpp"
 *     #undef main
 *
 * STRESS_SEED fixes the random seed, the default is fixed so CI runs are reproducible
 * STRESS_SCALE multiplies iteration counts for longer runs, the default is 1
 * The seed is printed to stderr so any failing run can be replayed
 *
***/

#include <bits/stdc++.h>

using namespace std;

namespace stress{
    inline long long env_or(const char* name, long long fallback){
        const char* value = getenv(name);
        return (value && *value) ? atoll(value) : fallback;
    }

    inline unsigned long long seed(){
        static const unsigned long long s = env_or("STRESS_SEED", 20260105);
        return s;
    }

    inline long long scaled(long long iterations){
        static const long long k = max(1LL, env_or("STRESS_SCALE", 1));
        return iterations * k;
    }

    inline mt19937_64& rng(){
        static mt19937_64 generator(seed());
        return generator;
    }

    inline long long rand_int(long long lo, long long hi){
        return uniform_int_distribution<long long>(lo, hi)(rng());
    }

    static const bool announced = fprintf(stderr, "STRESS_SEED=%llu\n", seed()) > 0;
}
