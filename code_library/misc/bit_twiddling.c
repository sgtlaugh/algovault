/// More bit twiddling hacks here - https://graphics.stanford.edu/~seander/bithacks.html

#include <stdio.h>
#include <assert.h>

/// Reverse the bits of an integer
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
unsigned int next_num(unsigned int x){
    assert(x > 0);
    unsigned int y = x & -x;
    x += y;
    unsigned int z = x & -x;
    z -= y;
    z = z >> bitscan(z & -z);
    return x | (z >> 1);
}

/// Returns the next number with the same number of 1 bits
unsigned int prev_num(unsigned int x){
    assert(x > 1);
    return ~next_num(~x);
}

/// Iterate over all the set bits of mask and store them in set_bits array
int iterate(unsigned long long mask, int* set_bits){
    int len = 0;
    while (mask){
        set_bits[len++] = __builtin_ctzll(mask);
        mask ^= (-mask & mask);
    }
    return len;
}

int main(){
    assert(reverse_bits(0) == 0);
    assert(reverse_bits(1) == 2147483648U);
    assert(reverse_bits(100) == 637534208);
    assert(reverse_bits(65536) == 32768);

    assert(bitscan(0) == 0);
    assert(bitscan(1) == 0);
    assert(bitscan(2) == 1);
    assert(bitscan(3) == 0);
    assert(bitscan(65536) == 16);

    assert(next_num(1) == 2);
    assert(next_num(16) == 32);
    assert(next_num(100) == 104);
    assert(next_num(2147483646) == 2684354559U);

    assert(prev_num(2) == 1);
    assert(prev_num(32) == 16);
    assert(prev_num(104) == 100);
    assert(prev_num(2684354559U) == 2147483646);

    int len, set_bits[64];

    len = iterate(0, set_bits);
    assert(len == 0);

    len = iterate(0b1000011010, set_bits);
    assert(len == 4);
    assert(set_bits[0] == 1 && set_bits[1] == 3 && set_bits[2] == 4 && set_bits[3] == 9);

    return 0;
}
