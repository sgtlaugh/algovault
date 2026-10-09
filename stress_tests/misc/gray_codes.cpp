#include "../common.h"

#define main library_main
#include "../../code_library/misc/gray_codes.c"
#undef main

/// Bit i of the gray code is bits i and i + 1 of x xored, bit i of the inverse is the xor of every bit from i upward
unsigned long long brute_gray(unsigned long long x){
    unsigned long long res = 0;
    for (int i = 0; i < 64; i++) res |= ((x >> i & 1) ^ (i < 63 ? x >> (i + 1) & 1 : 0)) << i;
    return res;
}

unsigned long long brute_inverse(unsigned long long g){
    unsigned long long res = 0, acc = 0;
    for (int i = 63; i >= 0; i--){
        acc ^= g >> i & 1;
        res |= acc << i;
    }
    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(500000); it++){
        unsigned long long x = stress::rng()() >> stress::rand_int(0, 63);
        if (it % 5 == 0) x = ~0ULL - stress::rand_int(0, 3);
        assert(gray_code(x) == brute_gray(x));
        assert(inverse_gray_code(x) == brute_inverse(x));
        assert(inverse_gray_code(gray_code(x)) == x && gray_code(inverse_gray_code(x)) == x);
        if (x != ~0ULL) assert(__builtin_popcountll(gray_code(x) ^ gray_code(x + 1)) == 1);  /// successive codes differ in one bit
    }

    return 0;
}
