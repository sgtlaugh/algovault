#include "../common.h"

#define main library_main
#include "../../code_library/misc/assembly.c"
#undef main

int random_int(){
    int x = (int)stress::rng()() >> stress::rand_int(0, 31);  /// both signs, small to full width
    int mode = stress::rand_int(0, 9);
    return mode == 0 ? INT_MIN : mode == 1 ? -1 : x;  /// INT_MIN / -1 is the case a signed division traps on
}

int main(){
    for (long long it = 0; it < stress::scaled(300000); it++){
        int x = (int)stress::rng()() >> stress::rand_int(0, 31);
        assert(popcount(x) == __builtin_popcount(x));
        assert(lzcount(x) == (x ? __builtin_clz(x) : 32));
        if (x) assert(bsf(x) == __builtin_ctz(x));

        int a = random_int(), b = random_int();
        if (it % 10 == 0) a = 0;
        assert(gcd(a, b) == (unsigned int)std::gcd(llabs(a), llabs(b)));  /// in long long, |INT_MIN| overflows int

        long double v = (long double)stress::rng()() / (stress::rand_int(1, 1000000) * 1.0L);
        assert(fsqrt(v) == sqrtl(v));
    }
    return 0;
}
