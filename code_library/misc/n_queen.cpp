/***
 *
 * N Queens Counting
 * Number of ways to place n non-attacking queens on an n x n board (https://oeis.org/A000170)
 *
 * Complexity: O(n!) worst case, bitmask backtracking over rows
 *
 * count_n_queens(n): number of solutions, 0 <= n <= 31 (asserted) because the column and diagonal masks are 32 bits
 * Tests check A000170 up to n = 15, a one-off run matched it up to n = 18
 * The running time grows about 7x per extra row, n = 17 takes about 30 seconds and n = 18 several minutes
 * The result is a long long because A000170 exceeds int from n = 19 on, that range is too slow to test
 *
 * Mirror symmetry: only first row columns in the left half are searched and doubled,
 * the middle column of an odd board is its own mirror and is counted once
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// unsigned: the diagonal masks shift past bit 30 from n = 16 on
long long n_queen_backtrack(int rows, unsigned int full, unsigned int cols, unsigned int left, unsigned int right){
    if (!rows) return 1;

    long long ways = 0;
    unsigned int free_cols = full & ~(cols | left | right);
    while (free_cols){
        unsigned int x = -free_cols & free_cols;
        free_cols ^= x;
        ways += n_queen_backtrack(rows - 1, full, cols | x, (left | x) << 1, (right | x) >> 1);
    }
    return ways;
}

long long count_n_queens(int n){
    assert(0 <= n && n <= 31);
    if (!n) return 1;

    unsigned int full = (1U << n) - 1;
    auto place = [&](int col){
        unsigned int x = 1U << col;
        return n_queen_backtrack(n - 1, full, x, x << 1, x >> 1);
    };

    long long ways = 0;
    for (int col = 0; col < n / 2; col++) ways += place(col);

    ways *= 2;
    if (n & 1) ways += place(n / 2);
    return ways;
}

int main(){
    assert(count_n_queens(0) == 1);
    assert(count_n_queens(1) == 1);
    assert(count_n_queens(2) == 0);
    assert(count_n_queens(3) == 0);
    assert(count_n_queens(4) == 2);
    assert(count_n_queens(5) == 10);
    assert(count_n_queens(6) == 4);
    assert(count_n_queens(8) == 92);
    assert(count_n_queens(13) == 73712);
    assert(count_n_queens(15) == 2279184);

    return 0;
}
