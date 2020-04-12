/***
 * iterates over all the set bits of a number, skipping 0 bits
 * also stores the set bits in an array passed as an argument
 * iterate2() should usually be preferred because iterate1() requires an expensive modulo operation
 * although the performance can vary and in some rare cases iterate() could be more desirable
 *
***/

#include <stdio.h>

#define MOD 67

char bithash[MOD];

/// call once for iterate1
void generate(){
    for (int i = 0; i < 64; i++) bithash[(1ULL << i) % MOD] = i;
}

/// iterate over all the set bits of mask and store them in set_bits array
int iterate1(unsigned long long mask, int* set_bits){
    int len = 0;
    while (mask){
        unsigned long long v = (mask & -mask);
        set_bits[len++] = bithash[v % MOD];
        mask ^= v;
    }
    return len;
}

/// iterate over all the set bits of mask and store them in set_bits array
int iterate2(unsigned long long mask, int* set_bits){
    int len = 0;
    while (mask){
        set_bits[len++] = __builtin_ctzll(mask);
        mask ^= (-mask & mask);
    }
    return len;
}

int main(){
    generate();
    int len, set_bits[64];
    unsigned long long x = 0b100000001011000100100101001010101001010110;

    len = iterate1(x, set_bits);
    printf("Set bits of %llu are:\n\n", x);
    for (int i = 0; i < len; i++){
        printf("Bit %d\n", set_bits[i]);
    }

    len = iterate2(x, set_bits);
    printf("\nSet bits of %llu are:\n\n", x);
    for (int i = 0; i < len; i++){
        printf("Bit %d\n", set_bits[i]);
    }

    return 0;
}
