#include "common.h"

#define main library_main
#include "../code_library/freivalds_algorithm.cpp"
#undef main

Matrix random_matrix(int row, int col, int range){
    Matrix res(row, col);
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++) res.mat[i][j] = stress::rand_int(-range, range);
    }
    return res;
}

int main(){
    static Matrix A, B, C, product;

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 8 : 30), k = stress::rand_int(1, it % 10 ? 8 : 30), m = stress::rand_int(1, it % 10 ? 8 : 30);
        int range = vector<int>{1, 1000, 500000000}[stress::rand_int(0, 2)];  /// the largest keeps every product entry within long long
        A = random_matrix(n, k, range), B = random_matrix(k, m, range);

        product = Matrix(n, m);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                for (int l = 0; l < k; l++) product.mat[i][j] += A.mat[i][l] * B.mat[l][j];
            }
        }
        assert(verify(A, B, product));

        /// A single wrong entry, off by as little as one
        C = product;
        C.mat[stress::rand_int(0, n - 1)][stress::rand_int(0, m - 1)] += stress::rand_int(0, 1) ? 1 : -stress::rand_int(1, 1000);
        assert(!verify(A, B, C));

        /// Errors that cancel in row sums, which an all-ones probe vector would miss
        if (m >= 2){
            C = product;
            int i = stress::rand_int(0, n - 1), j = stress::rand_int(0, m - 2);
            C.mat[i][j]++, C.mat[i][j + 1]--;
            assert(!verify(A, B, C));
        }

        /// Mismatched shapes
        C = Matrix(n, m + 1);
        assert(!verify(A, B, C));
        assert(!verify(A, Matrix(k + 1, m), product));
    }
    return 0;
}
