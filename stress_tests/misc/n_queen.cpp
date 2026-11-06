#include "../common.h"

#define main library_main
#include "../../code_library/misc/n_queen.cpp"
#undef main

/// Queens placed row by row on '.' squares with flag arrays for columns and both diagonals
long long brute(int row, const vector<string>& board, vector<bool>& col, vector<bool>& diag, vector<bool>& anti){
    int size = board.size();
    if (row == size) return 1;

    long long ways = 0;
    for (int c = 0; c < size; c++){
        if (board[row][c] != '.' || col[c] || diag[row + c] || anti[row - c + size]) continue;
        col[c] = diag[row + c] = anti[row - c + size] = true;
        ways += brute(row + 1, board, col, diag, anti);
        col[c] = diag[row + c] = anti[row - c + size] = false;
    }

    return ways;
}

long long brute(const vector<string>& board){
    int size = board.size();
    vector<bool> col(size), diag(2 * size), anti(2 * size + 1);
    return brute(0, board, col, diag, anti);
}

/// Every permutation of columns checked pairwise, independent of any row-by-row search
long long permutation_brute(const vector<string>& board){
    int size = board.size();
    vector<int> perm(size);
    iota(perm.begin(), perm.end(), 0);

    long long ways = 0;
    do {
        bool ok = true;
        for (int i = 0; i < size && ok; i++){
            if (board[i][perm[i]] != '.') ok = false;
            for (int j = 0; j < i && ok; j++){
                if (abs(perm[i] - perm[j]) == i - j) ok = false;
            }
        }
        ways += ok;
    } while (next_permutation(perm.begin(), perm.end()));

    return ways;
}

vector<string> random_board(int size, int blocked_percent){
    vector<string> board(size, string(size, '.'));
    for (auto& row : board){
        for (auto& cell : row){
            if (stress::rand_int(0, 99) < blocked_percent) cell = '*';
        }
    }

    return board;
}

/// A 31-queen solution (i, 2i + 1), (15 + i, 2i) for i < 15, then (30, 30), with a few random extra '.' squares
vector<string> sparse_board_31(){
    vector<int> queen(31);
    for (int i = 0; i < 15; i++) queen[i] = 2 * i + 1, queen[15 + i] = 2 * i;
    queen[30] = 30;

    vector<string> board(31, string(31, '*'));
    for (int i = 0; i < 31; i++){
        board[i][queen[i]] = '.';
        if (stress::rand_int(0, 3) == 0) board[i][stress::rand_int(0, 30)] = '.';
    }

    return board;
}

/// count_n_queens(n) against A000170 and the brute force, count_n_queens(board) against both brute forces
int main(){
    const long long oeis[] = {1, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724, 2680, 14200, 73712, 365596};  /// A000170
    for (int size = 0; size <= 14; size++){
        if (size <= 10) assert(brute(vector<string>(size, string(size, '.'))) == oeis[size]);
        assert(count_n_queens(size) == oeis[size]);
        assert(count_n_queens(vector<string>(size, string(size, '.'))) == oeis[size]);
    }

    /// Every blocked pattern of the boards up to 4 x 4
    for (int size = 1; size <= 4; size++){
        for (int mask = 0; mask < (1 << (size * size)); mask++){
            vector<string> board(size, string(size, '.'));
            for (int cell = 0; cell < size * size; cell++){
                if (mask >> cell & 1) board[cell / size][cell % size] = '*';
            }
            assert(count_n_queens(board) == permutation_brute(board));
        }
    }

    /// Random blocked squares, any non-'.' character blocks
    for (int it = 0; it < stress::scaled(400); it++){
        int size = stress::rand_int(5, 8);
        vector<string> board = random_board(size, stress::rand_int(0, 40));
        if (stress::rand_int(0, 1)) board[stress::rand_int(0, size - 1)][stress::rand_int(0, size - 1)] = 'Q';
        assert(count_n_queens(board) == permutation_brute(board));
    }

    for (int it = 0; it < stress::scaled(40); it++){
        int size = stress::rand_int(9, 12);
        vector<string> board = random_board(size, stress::rand_int(0, 30));
        assert(count_n_queens(board) == brute(board));
    }

    /// n = 31 uses bit 31 of every mask, a signed shift or a 31-bit full mask breaks it
    for (int it = 0; it < stress::scaled(40); it++){
        vector<string> board = sparse_board_31();
        long long ways = count_n_queens(board);
        assert(ways >= 1);
        assert(ways == brute(board));
    }

    stress::expect_assert("misc/n_queen.cpp", []{ count_n_queens(32); });
    stress::expect_assert("misc/n_queen.cpp", []{ count_n_queens(vector<string>(32, string(32, '.'))); });
    stress::expect_assert("misc/n_queen.cpp", []{ count_n_queens(vector<string>{"...", "..", "..."}); });

    return 0;
}
