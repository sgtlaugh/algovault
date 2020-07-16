#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAX 101

using namespace std;

struct Matrix{
    int row, col;
    long long mat[MAX][MAX];

    Matrix(){}
    Matrix(int n, int m, int diagonal = 0){
        memset(mat, 0, sizeof(mat));
        row = n, col = m;
        for (int i = min(n, m) - 1; i >= 0; i--) mat[i][i] = diagonal;
    }

    Matrix operator* (const Matrix& other) const{
        int i, j, k;
        Matrix res(row, other.col);

        for(i = 0; i < row; i++){
            for(j = 0; j < other.col; j++){
                res.mat[i][j] = 0;
                for(k = 0; k < col; k++){
                    res.mat[i][j] += mat[i][k] * other.mat[k][j];
                }
            }
        }

        return res;
    }

    bool operator == (const Matrix& other) const{
        if (row != other.row || col != other.col) return false;

        for (int i = 0; i < row; i++){
            for (int j = 0; j < col; j++){
                if (mat[i][j] != other.mat[i][j]) return false;
            }
        }
        return true;
    }

    bool operator != (const Matrix & other) const {
        return !(*this == other);
    }
};

/***
 *
 * a randomized algorithm that verifies matrix multiplication
 * returns if A * B = C with probability of failure less than 1/(2^ntrials)
 * complexity: O(n^3) * ntrials
 *
***/

bool verify(Matrix A, Matrix B, Matrix C, int ntrials=30){
    if (A.col != B.row || C.row != A.row || C.col != B.col) return false;

    for (int l = 0; l < ntrials; l++){
        Matrix D = Matrix(C.col, 1);
        for (int i = 0; i < D.row; i++) D.mat[i][0] = rand();

        Matrix X = A * (B * D);
        Matrix Y = C * D;
        if (X != Y) return false;
    }

    return true;
}

int main(){
    Matrix A = Matrix(2, 3);
    Matrix B = Matrix(3, 4);

    for (int i = 0; i < A.row; i++){
        for (int j = 0; j < A.col; j++){
            A.mat[i][j] = rand() % 6 - 3;
        }
    }

    for (int i = 0; i < B.row; i++){
        for (int j = 0; j < B.col; j++){
            B.mat[i][j] = rand() % 6 - 3;
        }
    }

    Matrix C = A * B;
    assert(verify(A, B, C));

    for (int i = 0; i < C.row; i++){
        for (int j = 0; j < C.col; j++){
            C.mat[i][j]--;
            assert(!verify(A, B, C));
            C.mat[i][j] += 2;
            assert(!verify(A, B, C));
            C.mat[i][j]--;
        }
    }

    return 0;
}
