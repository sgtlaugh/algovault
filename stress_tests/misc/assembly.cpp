#include "../common.h"

#define main library_main
#include "../../code_library/misc/assembly.c"
#undef main

int random_int(){
    int x = stress::rng()() >> 33;  /// non-negative, the gcd's division zero-extends its dividend
    return x >> stress::rand_int(0, 30);
}

int main(){
    for (long long it = 0; it < stress::scaled(300000); it++){
        int x = (int)stress::rng()() >> stress::rand_int(0, 31);
        assert(popcount(x) == __builtin_popcount(x));
        assert(lzcount(x) == (x ? __builtin_clz(x) : 32));
        if (x) assert(bsf(x) == __builtin_ctz(x));

        int a = random_int(), b = random_int();
        if (it % 10 == 0) a = 0;
        assert(gcd(a, b) == std::gcd(a, b));

        long double v = (long double)stress::rng()() / (stress::rand_int(1, 1000000) * 1.0L);
        assert(fsqrt(v) == sqrtl(v));
    }
    return 0;
}
