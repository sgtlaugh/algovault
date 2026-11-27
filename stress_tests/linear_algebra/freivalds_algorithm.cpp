#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/freivalds_algorithm.cpp"
#undef main

FreivaldsMatrix random_matrix(int row, int col, int range){
    FreivaldsMatrix res(row, col);
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++) res[i][j] = stress::rand_int(-range, range);
    }

    return res;
}

int main(){
    static FreivaldsMatrix A, B, C, product;

    for (long long it = 0; it < stress::scaled(1500); it++){
        int lim = it % 100 == 0 ? 150 : it % 10 ? 8 : 30;  /// 150 is past any fixed 101 x 101 storage
        int n = stress::rand_int(1, lim), k = stress::rand_int(1, lim), m = stress::rand_int(1, lim);
        /// The largest range keeps every product entry within long long
        int range = vector<int>{1, 1000, lim > 30 ? 200000000 : 500000000}[stress::rand_int(0, 2)];
        A = random_matrix(n, k, range), B = random_matrix(k, m, range);

        product = FreivaldsMatrix(n, m);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                for (int l = 0; l < k; l++) product[i][j] += A[i][l] * B[l][j];
            }
        }
        assert(freivalds_verify(A, B, product) && A * B == product);

        /// A single wrong entry, off by as little as one
        C = product;
        C[stress::rand_int(0, n - 1)][stress::rand_int(0, m - 1)] += stress::rand_int(0, 1) ? 1 : -stress::rand_int(1, 1000);
        assert(!freivalds_verify(A, B, C));

        /// Errors that cancel in row sums, which an all-ones probe vector would miss
        if (m >= 2){
            C = product;
            int i = stress::rand_int(0, n - 1), j = stress::rand_int(0, m - 2);
            C[i][j]++, C[i][j + 1]--;
            assert(!freivalds_verify(A, B, C));
        }

        /// Mismatched shapes
        C = FreivaldsMatrix(n, m + 1);
        assert(!freivalds_verify(A, B, C));
        assert(!freivalds_verify(A, FreivaldsMatrix(k + 1, m), product));
    }

    /// A row of ones times (LLONG_MAX, 1, -1) is LLONG_MAX: its partial sums overflow, the entry does not
    FreivaldsMatrix ones(1, 3), column(3, 1);
    for (int k = 0; k < 3; k++) ones[0][k] = 1;
    column[0][0] = LLONG_MAX, column[1][0] = 1, column[2][0] = -1;
    assert((ones * column)[0][0] == LLONG_MAX);

    /// A long long diagonal is kept whole, an int parameter cut 5e9 to 705032704
    long long diagonal = stress::rand_int(5000000000LL, 5000000000LL);
    FreivaldsMatrix big(2, 2, diagonal);
    assert(big[0][0] == 5000000000LL && big[1][1] == 5000000000LL && big[0][1] == 0);

    return 0;
}
