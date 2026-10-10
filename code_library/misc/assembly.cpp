/***
 *
 * Inline x86 Assembly
 * Examples of GCC extended asm: popcount, leading zero count, bit scan, gcd and x87 square root
 *
 * Complexity: O(1) per call, O(log min(|a|, |b|)) for asm_gcd
 *
 * Requires x86-64 and GCC or Clang (AT&T syntax)
 * POPCNT and LZCNT need CPU support, without LZCNT the instruction silently runs as BSR and returns a wrong count
 * For real use prefer __builtin_popcount, __builtin_clz, __builtin_ctz and std::gcd: portable, and std::gcd is ~1.6x faster
 *
 * asm_gcd works on magnitudes and returns unsigned, so asm_gcd(INT_MIN, 0) = 2^31
 * asm_bsf(0) returns 0, asm_lzcnt(0) returns 32
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

int asm_bsf(unsigned int x){
    if (!x) return 0;  /// bsf leaves its destination undefined for 0

    __asm__("bsf %0, %0" : "=r" (x) : "0" (x));
    return x;
}

/// divl on magnitudes: idivl traps on INT_MIN / -1, and gcd(INT_MIN, 0) = 2^31 needs unsigned
unsigned int asm_gcd(int a, int b){
    unsigned int x = a < 0 ? 0u - a : (unsigned int)a, y = b < 0 ? 0u - b : (unsigned int)b, res;

    __asm__(
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

int asm_lzcnt(int x){
    int counter = 0;
    __asm__("LZCNT %1, %0;" : "=r"(counter) : "r"(x));
    return counter;
}

int asm_popcount(int x){
    int counter = 0;
    __asm__("POPCNT %1, %0;" : "=r"(counter) : "r"(x));
    return counter;
}

long double asm_sqrt(long double x){
    __asm__("fsqrt" : "+t" (x));
    return x;
}

int main(){
    assert(asm_popcount(13) == 3);               /// 1101
    assert(asm_popcount(-1) == 32);
    assert(asm_lzcnt(100) == 25);                /// 100 needs 7 bits
    assert(asm_lzcnt(0) == 32);
    assert(asm_bsf(100) == 2);                   /// 1100100
    assert(asm_gcd(1071, 462) == 21);
    assert(asm_gcd(-12, -18) == 6);              /// signs are dropped
    assert(asm_gcd(INT_MIN, 0) == 2147483648u);  /// the reason the result is unsigned
    assert(asm_sqrt(16.0L) == 4.0L);
    return 0;
}
