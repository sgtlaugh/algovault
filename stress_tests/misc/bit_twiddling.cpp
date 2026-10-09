#include "../common.h"

#define main library_main
#include "../../code_library/misc/bit_twiddling.cpp"
#undef main

/// Gosper's hack, an independent formula for the next number with the same popcount
unsigned int gosper(unsigned int x){
    unsigned int c = x & -x, r = x + c;
    return (((r ^ x) >> 2) / c) | r;
}

bool has_next(unsigned int x){
    unsigned int top = __builtin_popcount(x) ? ~0U << (32 - __builtin_popcount(x)) : 0;
    return x && x != top;  /// all set bits packed at the top have no successor
}

bool has_prev(unsigned int x){
    return x != ~0U && (x & (x + 1));  /// all set bits packed at the bottom have no predecessor
}

unsigned int random_word(){
    unsigned int x = stress::rng()();
    int mode = stress::rand_int(0, 3);
    if (mode == 0) x >>= stress::rand_int(0, 31);
    if (mode == 1) x &= stress::rng()() & stress::rng()();  /// sparse
    if (mode == 2) x |= stress::rng()() | stress::rng()();  /// dense
    return x;
}

void check_neighbours(unsigned int x){
    if (has_next(x)){
        unsigned int y = next_same_popcount(x);
        assert(y == gosper(x) && prev_same_popcount(y) == x);
    }
    if (has_prev(x)){
        unsigned int y = prev_same_popcount(x);
        assert(y < x && __builtin_popcount(y) == __builtin_popcount(x) && next_same_popcount(y) == x);
    }
}

int main(){
    /// Brute force neighbours on small values
    for (unsigned int x = 1; x < (1u << 15); x++){
        unsigned int y = x + 1;
        while (__builtin_popcount(y) != __builtin_popcount(x)) y++;
        assert(next_same_popcount(x) == y);
        if (x >= 2 && has_prev(x)){
            unsigned int z = x - 1;
            while (__builtin_popcount(z) != __builtin_popcount(x)) z--;
            assert(prev_same_popcount(x) == z);
        }
    }

    /// Every word with popcount 1, 2, 30 or 31, which includes all the boundary words without a neighbour
    for (int i = 0; i < 32; i++){
        for (int j = i; j < 32; j++){
            unsigned int w = (1u << i) | (1u << j);
            check_neighbours(w);
            check_neighbours(~w);
        }
    }

    for (long long it = 0; it < stress::scaled(300000); it++){
        unsigned int x = random_word();

        unsigned int r = 0;
        for (int i = 0; i < 32; i++) r |= (x >> i & 1) << (31 - i);
        assert(reverse_bits(x) == r);

        check_neighbours(x);

        unsigned long long mask = (unsigned long long)random_word() << stress::rand_int(0, 32) | random_word();
        vector<int> bits = set_bit_indices(mask);
        assert((int)bits.size() == __builtin_popcountll(mask));
        for (int i = 0, j = 0; i < 64; i++){
            if (mask >> i & 1) assert(bits[j++] == i);
        }
    }

    return 0;
}
