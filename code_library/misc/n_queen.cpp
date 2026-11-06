/***
 *
 * N Queens Counting
 * Number of ways to place n non-attacking queens on an n x n board, optionally with blocked squares
 *
 * Complexity: O(n!) worst case, bitmask backtracking over rows
 *
 * count_n_queens(n): number of solutions on the empty board (https://oeis.org/A000170), 0 <= n <= 31
 * count_n_queens(board): number of solutions with queens only on '.' squares of n strings of length n, n <= 31
 *
 * Notes:
 *   - Any character other than '.' is a blocked square, which holds no queen but does not stop an attack through it
 *   - The running time grows about 7x per extra row, n = 17 takes about 30 seconds and n = 18 several minutes
 *   - count_n_queens(board) is about 2x slower than count_n_queens(n) on an empty board, it cannot use mirror symmetry
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// unsigned: the diagonal masks shift past bit 30 from n = 16 on
template <typename Usable>  /// usable(rows) is a functor so count_n_queens(n) keeps its mask in a register, 8% faster
long long n_queen_backtrack(const Usable& usable, int rows, unsigned int cols, unsigned int left, unsigned int right){
    if (!rows) return 1;

    long long ways = 0;
    unsigned int free_cols = usable(rows) & ~(cols | left | right);
    while (free_cols){
        unsigned int x = -free_cols & free_cols;
        free_cols ^= x;
        ways += n_queen_backtrack(usable, rows - 1, cols | x, (left | x) << 1, (right | x) >> 1);
    }

    return ways;
}

long long count_n_queens(int n){
    assert(0 <= n && n <= 31);
    if (!n) return 1;

    unsigned int full = (1U << n) - 1;
    auto usable = [full](int){ return full; };
    auto place = [&](int col){
        unsigned int x = 1U << col;
        return n_queen_backtrack(usable, n - 1, x, x << 1, x >> 1);
    };

    long long ways = 0;
    for (int col = 0; col < n / 2; col++) ways += place(col);

    ways *= 2;
    if (n & 1) ways += place(n / 2);
    return ways;
}

long long count_n_queens(const vector<string>& board){
    int n = board.size();
    assert(n <= 31);

    vector<unsigned int> usable(n);
    for (int i = 0; i < n; i++){
        assert((int)board[i].size() == n);
        for (int j = 0; j < n; j++){
            if (board[i][j] == '.') usable[i] |= 1U << j;
        }
    }

    return n_queen_backtrack([&](int rows){ return usable[n - rows]; }, n, 0, 0, 0);
}

int main(){
    /// The classic 8 queens puzzle has 92 solutions, the 3 x 3 board none
    assert(count_n_queens(8) == 92);
    assert(count_n_queens(4) == 2);
    assert(count_n_queens(3) == 0);

    /// The two 4 x 4 solutions are .Q.. / ...Q / Q... / ..Q. and its mirror ..Q. / Q... / ...Q / .Q..
    assert(count_n_queens({".*..", "....", "....", "...."}) == 1);
    assert(count_n_queens({"....", ".*..", "....", "...."}) == 2);

    return 0;
}
