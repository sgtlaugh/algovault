#include "../common.h"

#define main library_main
#include "../../code_library/misc/n_queen.cpp"
#undef main

/// Queens placed row by row with flag arrays for columns and both diagonals
long long brute(int row, int size, vector<bool>& col, vector<bool>& diag, vector<bool>& anti){
    if (row == size) return 1;

    long long ways = 0;
    for (int c = 0; c < size; c++){
        if (col[c] || diag[row + c] || anti[row - c + size]) continue;
        col[c] = diag[row + c] = anti[row - c + size] = true;
        ways += brute(row + 1, size, col, diag, anti);
        col[c] = diag[row + c] = anti[row - c + size] = false;
    }
    return ways;
}

long long brute(int size){
    vector<bool> col(size), diag(2 * size), anti(2 * size + 1);
    return brute(0, size, col, diag, anti);
}

/// count_n_queens against the row-by-row brute force for n <= 10 and against A000170 up to n = 14
int main(){
    const long long oeis[] = {1, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724, 2680, 14200, 73712, 365596};  /// A000170
    for (int size = 0; size <= 14; size++){
        if (size <= 10) assert(brute(size) == oeis[size]);
        assert(count_n_queens(size) == oeis[size]);
    }

    return 0;
}
