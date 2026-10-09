/***
 * 
 * https://en.wikipedia.org/wiki/Eight_queens_puzzle
 * 
 * Counts the number of ways for the N queens puzzle
 * Mirror symmetry: only first row columns in the left half are searched and doubled,
 * the middle column of an odd board is its own mirror and is counted once
 * 
***/

#include <stdio.h>
#include <assert.h>

int n, counter;
unsigned int lim;

void backtrack(int i, unsigned int c, unsigned int l, unsigned int r){  /// unsigned, diagonal bits shift past bit 30 from n = 16 on
    if (!i){
        counter++;
        return;
    }

    unsigned int bitmask, x;
    --i, bitmask = lim & ~(l | r | c);

    while (bitmask){
        x = (-bitmask & bitmask);
        bitmask ^= x;
        backtrack(i, c | x, (l | x) << 1, (r | x) >> 1);
    }
}

int count_ways(int dimension){
    n = dimension;
    if (!n) return 1;

    int i;
    counter = 0, lim = (1U << n) - 1;
    for (i = 0; i < n / 2; i++) backtrack(n - 1, 1U << i, 2U << i, (1U << i) >> 1);

    counter *= 2;
    if (n & 1) backtrack(n - 1, 1U << i, 2U << i, (1U << i) >> 1);
    return counter;
}

int main(){
    assert(count_ways(1) == 1);
    assert(count_ways(2) == 0);
    assert(count_ways(3) == 0);
    assert(count_ways(4) == 2);
    assert(count_ways(8) == 92);
    assert(count_ways(13) == 73712);
    assert(count_ways(15) == 2279184);

    return 0;
}
