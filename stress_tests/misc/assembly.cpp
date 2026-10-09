#include "../common.h"

#define main library_main
#include "../../code_library/misc/assembly.cpp"
#undef main

int random_int(){
    int x = (int)stress::rng()() >> stress::rand_int(0, 31);  /// both signs, small to full width
    int mode = stress::rand_int(0, 9);
    return mode == 0 ? INT_MIN : mode == 1 ? -1 : x;  /// INT_MIN / -1 is the case a signed division traps on
}

unsigned int reference_gcd(int a, int b){
    return (unsigned int)std::gcd(llabs(a), llabs(b));  /// in long long, |INT_MIN| overflows int
}

/// compares every asm routine against the GCC builtins, std::gcd and sqrtl
int main(){
    for (int a = -60; a <= 60; a++){
        for (int b = -60; b <= 60; b++) assert(asm_gcd(a, b) == reference_gcd(a, b));
    }

    for (int x = -70000; x <= 70000; x++){
        assert(asm_popcount(x) == __builtin_popcount(x));
        assert(asm_lzcnt(x) == (x ? __builtin_clz(x) : 32));
        assert(asm_bsf(x) == (x ? __builtin_ctz(x) : 0));
    }

    for (long long r = 0; r <= 100000; r++) assert(asm_sqrt((long double)(r * r)) == (long double)r);

    for (long long it = 0; it < stress::scaled(300000); it++){
        int x = (int)stress::rng()() >> stress::rand_int(0, 31);
        assert(asm_popcount(x) == __builtin_popcount(x));
        assert(asm_lzcnt(x) == (x ? __builtin_clz(x) : 32));
        assert(asm_bsf(x) == (x ? __builtin_ctz(x) : 0));

        int a = random_int(), b = random_int();
        if (it % 10 == 0) a = 0;
        assert(asm_gcd(a, b) == reference_gcd(a, b));

        long double v = (long double)stress::rng()() / (stress::rand_int(1, 1000000) * 1.0L);
        assert(asm_sqrt(v) == sqrtl(v));
    }

    return 0;
}
