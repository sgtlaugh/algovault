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
    int t, n;
    unsigned long long k;

    scanf("%d", &t);
    while (t--){
        scanf("%d %llu", &n, &k);
        printf("%llu\n", gray_code(k));
    }

    return 0;
}
