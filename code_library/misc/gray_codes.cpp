/***
 *
 * Gray Codes
 * Binary reflected Gray code of a 64-bit unsigned integer and its inverse
 *
 * Complexity: O(1) for gray_code, O(log w) = 6 shifts for inverse_gray_code
 *
 * Successive codes differ in exactly one bit: the 3-bit sequence is 000, 001, 011, 010, 110, 111, 101, 100
 * inverse_gray_code(gray_code(x)) == x for every 64-bit x
 * https://en.wikipedia.org/wiki/Gray_code
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

unsigned long long gray_code(unsigned long long x){
    return x ^ (x >> 1);
}

/// Prefix xor of all higher bits by doubling, 6 steps instead of one per bit
unsigned long long inverse_gray_code(unsigned long long x){
    x ^= x >> 1;
    x ^= x >> 2;
    x ^= x >> 4;
    x ^= x >> 8;
    x ^= x >> 16;
    x ^= x >> 32;
    return x;
}

int main(){
    for (unsigned long long x = 0; x < 1048576; x++){
        assert(inverse_gray_code(gray_code(x)) == x);
    }

    vector<unsigned long long> three_bit = {0, 1, 3, 2, 6, 7, 5, 4};
    for (unsigned long long x = 0; x < 8; x++){
        assert(gray_code(x) == three_bit[x]);
        assert(inverse_gray_code(three_bit[x]) == x);
    }

    assert(gray_code(1000000007) == 643280644);
    assert(gray_code(1000000000000000003ULL) == 797398725282889730ULL);
    assert(gray_code(~0ULL) == 1ULL << 63);

    assert(inverse_gray_code(643280644) == 1000000007);
    assert(inverse_gray_code(797398725282889730ULL) == 1000000000000000003ULL);
    assert(inverse_gray_code(1ULL << 63) == ~0ULL);

    return 0;
}
