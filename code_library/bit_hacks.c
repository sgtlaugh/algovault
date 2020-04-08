#include <stdio.h>
#include <assert.h>


/// More bit twiddling hacks here - https://graphics.stanford.edu/~seander/bithacks.html

unsigned int reverse_bits(unsigned int v){
    v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    return ((v >> 16) | (v << 16));
}

/// Returns i if x = 2^i and 0 otherwise
unsigned int bitscan(unsigned int x){
    __asm__ volatile("bsf %0, %0" : "=r" (x) : "0" (x));
    return x;
}

/// Returns the next number with the same number of 1 bits
unsigned int nextnum(unsigned int x){
    assert(x > 0);
    unsigned int y = x & -x;
    x += y;
    unsigned int z = x & -x;
    z -= y;
    z = z >> bitscan(z & -z);
    return x | (z >> 1);
}

/// Returns the next number with the same number of 1 bits
unsigned int prevnum(unsigned int x){
    assert(x > 1);
    return ~nextnum(~x);
}

int main(){
    printf("%u\n", reverse_bits(0));      /// 0
    printf("%u\n", reverse_bits(1));      /// 2147483648
    printf("%u\n", reverse_bits(100));    /// 637534208
    printf("%u\n", reverse_bits(65536));  /// 32768

    printf("%u\n", bitscan(0));           /// 0
    printf("%u\n", bitscan(1));           /// 0
    printf("%u\n", bitscan(2));           /// 1
    printf("%u\n", bitscan(3));           /// 0
    printf("%u\n", bitscan(65536));       /// 16

    printf("%u\n", nextnum(1));           /// 2
    printf("%u\n", nextnum(16));          /// 32
    printf("%u\n", nextnum(100));         /// 104
    printf("%u\n", nextnum(2147483646));  /// 2684354559

    printf("%u\n", prevnum(2));           /// 1
    printf("%u\n", prevnum(32));          /// 16
    printf("%u\n", prevnum(104));         /// 100
    printf("%u\n", prevnum(2684354559U)); /// 2147483646

    return 0;
}
