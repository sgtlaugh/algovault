#include "../common.h"

#define main library_main
#include "../../code_library/misc/n_queen.c"
#undef main

/// Queens placed row by row with flag arrays for columns and both diagonals
int brute(int row, int size, vector<bool>& col, vector<bool>& diag, vector<bool>& anti){
    if (row == size) return 1;
    int ways = 0;
    for (int c = 0; c < size; c++){
        if (col[c] || diag[row + c] || anti[row - c + size]) continue;
        col[c] = diag[row + c] = anti[row - c + size] = true;
        ways += brute(row + 1, size, col, diag, anti);
        col[c] = diag[row + c] = anti[row - c + size] = false;
    }
    return ways;
}

int main(){
    const int oeis[] = {1, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724, 2680, 14200};  /// A000170
    for (long long it = 0; it < stress::scaled(30); it++){
        int size = stress::rand_int(0, it % 3 ? 9 : 12);
        int expected = oeis[size];
        if (size <= 10){
            vector<bool> col(size), diag(2 * size), anti(2 * size + 1);
            expected = brute(0, size, col, diag, anti);
            assert(expected == oeis[size]);
        }
        assert(count_ways(size) == expected);
    }
    return 0;
}
