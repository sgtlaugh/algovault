/***
 *
 * Inline x86 assembly examples: popcount, leading zero count, bit scan, gcd and x87 square root
 *
 * POPCNT and LZCNT need CPU support, without LZCNT the instruction silently runs as BSR and returns a wrong count
 * For real use prefer __builtin_popcount, __builtin_clz, __builtin_ctz and std::gcd: portable, and std::gcd is ~1.6x faster
 *
***/

#include <stdio.h>
#include <assert.h>

int popcount(int x){
    int counter = 0;
    __asm__ volatile("POPCNT %1, %0;"
        :"=r"(counter)
        :"r"(x)
        :
    );
    return counter;
}

int lzcount(int x){
    int counter = 0;
    __asm__ volatile("LZCNT %1, %0;"
        :"=r"(counter)
        :"r"(x)
        :
    );
    return counter;
}

int bsf(unsigned int x){
    __asm__ volatile("bsf %0, %0" : "=r" (x) : "0" (x));
    return x;
}

/// divl on magnitudes: idivl traps on INT_MIN / -1, and gcd(INT_MIN, 0) = 2^31 needs unsigned
unsigned int gcd(int a, int b){
    unsigned int x = a < 0 ? 0u - a : (unsigned int)a, y = b < 0 ? 0u - b : (unsigned int)b, res;
    __asm__ volatile(
                      "movl %1, %%eax;"
                      "movl %2, %%ebx;"
                      "repeat_%=:\n"
                      "cmpl $0, %%ebx;"
                      "je terminate_%=\n;"
                      "xorl %%edx, %%edx;"
                      "divl %%ebx;"
                      "movl %%ebx, %%eax;"
                      "movl %%edx, %%ebx;"
                      "jmp repeat_%=\n;"
                      "terminate_%=:\n"
                      "movl %%eax, %0;"

                      : "=g"(res)
                      : "g"(x), "g"(y)
                      : "eax", "ebx", "edx"

    );

    return res;
}

long double fsqrt(long double x){
    __asm__ volatile("fsqrt" : "+t" (x));
    return x;
}

int main(){
    assert(popcount(0) == 0);
    assert(popcount(1) == 1);
    assert(popcount(13) == 3);
    assert(popcount(2147483647) == 31);

    assert(lzcount(0) == 32);
    assert(lzcount(100) == 25);

    assert(bsf(0) == 0);
    assert(bsf(15) == 0);
    assert(bsf(16) == 4);
    assert(bsf(100) == 2);

    assert(gcd(0, 10) == 10);
    assert(gcd(120, 1260) == 60);
    assert(gcd(-4, 6) == 2 && gcd(-4, 1) == 1 && gcd(-12, -18) == 6);
    assert(gcd(-2147483647 - 1, -1) == 1 && gcd(-2147483647 - 1, 0) == 2147483648u);

    return 0;
}
