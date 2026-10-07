/***
 * 
 * https://en.wikipedia.org/wiki/Eight_queens_puzzle
 * 
 * Counts the number of ways for the N queens puzzle
 * Can be further optimized if we consider board symmetry
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
    counter = 0, lim = (1U << n) - 1;
    backtrack(n, 0, 0, 0);

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
