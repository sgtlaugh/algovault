#include <stdio.h>

#define MOD 67


char bithash[MOD];
int len, set_bits[64];

/// call once
void generate(){
    int i;
    for (i = 0; i < 64; i++) bithash[(1ULL << i) % MOD] = i;
}

/// iterate over all the set bits of mask and store them in set_bits array
void iterate(unsigned long long mask){
    len = 0;

    while (mask){
        unsigned long long v = (mask & -mask);
        set_bits[len++] = bithash[v % MOD];
        mask ^= v;
    }
}

int main(){
    generate();

    unsigned long long x = 0b100000001011000100100101001010101001010110;
    iterate(x);

    int i;
    printf("Set bits of %llu are:\n\n", x);
    for (i = 0; i < len; i++){
        printf("Bit %d\n", set_bits[i]);
    }

    return 0;
}
