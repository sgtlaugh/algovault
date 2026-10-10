/***
 *
 * Bit Twiddling
 * Reverse bits, step to the next or previous word with the same popcount, list set bits
 *
 * Complexity: O(1) per operation, O(popcount) for set_bit_indices
 *
 * reverse_bits(x)            bit i of x moves to bit 31 - i
 * next_same_popcount(x)      smallest y > x with popcount(y) == popcount(x)
 *                            requires a successor: x != 0 and the set bits of x not all packed at the top
 * prev_same_popcount(x)      largest y < x with popcount(y) == popcount(x)
 *                            requires a predecessor: x != 0 and the set bits of x not all packed at the bottom
 * set_bit_indices(mask)      indices of the set bits of a 64-bit mask, ascending
 *
 * The index of the lowest set bit is __builtin_ctz(x) (undefined for x = 0)
 * More tricks: https://graphics.stanford.edu/~seander/bithacks.html
 *
 * Example:
 *   for (unsigned x = (1u << k) - 1; x < (1u << n); x = next_same_popcount(x))  // k-subsets, 1 <= k <= n < 32
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

unsigned int reverse_bits(unsigned int v){
    v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    return (v >> 16) | (v << 16);
}

unsigned int next_same_popcount(unsigned int x){
    unsigned int y = x & -x;
    assert(x + y != 0);  /// x == 0 or bits packed at the top: no successor, x + y wraps to 0

    x += y;
    unsigned int z = (x & -x) - y;
    z >>= __builtin_ctz(z);

    return x | (z >> 1);
}

unsigned int prev_same_popcount(unsigned int x){
    return ~next_same_popcount(~x);
}

vector<int> set_bit_indices(unsigned long long mask){
    vector<int> indices;
    for (; mask; mask &= mask - 1) indices.push_back(__builtin_ctzll(mask));
    return indices;
}

int main(){
    assert(reverse_bits(1) == 0x80000000U);
    assert(reverse_bits(0x0000FFFFU) == 0xFFFF0000U);

    /// The next and previous word with the same three bits set: 0111 -> 1011 -> 1101
    assert(next_same_popcount(0b0111) == 0b1011);
    assert(next_same_popcount(0b1011) == 0b1101);
    assert(prev_same_popcount(0b1011) == 0b0111);

    assert((set_bit_indices(0b1000011010) == vector<int>{1, 3, 4, 9}));
    assert((set_bit_indices(1ULL << 63) == vector<int>{63}));
    return 0;
}
