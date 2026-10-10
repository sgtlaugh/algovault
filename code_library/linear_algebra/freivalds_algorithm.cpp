/***
 *
 * A randomized algorithm that verifies matrix multiplication
 * Checks and returns if A * B = C with probability of failure less than 1 / (2^number_of_trials)
 * The check runs modulo 2^64, which is exact whenever C and the true product A * B fit in long long
 * Complexity: O(n^2) * number_of_trials, where n is the size of the matrix, O(row * col) memory per Matrix
 *
 * Matrix M(row, col, diagonal = 0) is a row x col matrix of long long, M[i][j] for 0 <= i < row, 0 <= j < col
 * verify(A, B, C, ntrials = 30) returns false on mismatched shapes
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/// One row-major vector: vector<vector<long long>> measured 7-11% slower on a 100 x 100 product
struct Matrix{
    int row, col;
    vector<long long> mat;

    Matrix() : row(0), col(0) {}
    Matrix(int row, int col, int diagonal = 0) : row(row), col(col), mat((size_t)row * col, 0){
        for (int i = min(row, col) - 1; i >= 0; i--) (*this)[i][i] = diagonal;
    }

    long long* operator[](int i){
        return mat.data() + (size_t)i * col;
    }

    const long long* operator[](int i) const{
        return mat.data() + (size_t)i * col;
    }

    Matrix operator* (const Matrix& other) const{
        Matrix res(row, other.col);
        for (int i = 0; i < row; i++){
            for (int j = 0; j < other.col; j++){
                long long sum = 0;
                for (int k = 0; k < col; k++) sum += (*this)[i][k] * other[k][j];
                res[i][j] = sum;
            }
        }

        return res;
    }

    bool operator == (const Matrix& m) const{
        return row == m.row && col == m.col && mat == m.mat;
    }

    bool operator != (const Matrix & m) const {
        return !(*this == m);
    }
};

unsigned long long dot(const long long* a, const vector<unsigned long long>& b){
    unsigned long long res = 0;
    for (size_t i = 0; i < b.size(); i++) res += (unsigned long long)a[i] * b[i];
    return res;
}

bool verify(const Matrix& A, const Matrix& B, const Matrix& C, int ntrials=30){
    if (A.col != B.row || C.row != A.row || C.col != B.col) return false;

    static mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    vector<unsigned long long> r(C.col), br(B.row), x(A.row), y(C.row);  /// unsigned: wraps mod 2^64 instead of overflowing

    for (int l = 0; l < ntrials; l++){
        for (auto& v : r) v = rng();
        for (int i = 0; i < B.row; i++) br[i] = dot(B[i], r);
        for (int i = 0; i < A.row; i++) x[i] = dot(A[i], br), y[i] = dot(C[i], r);
        if (x != y) return false;
    }

    return true;
}

int main(){
    srand(0);
    Matrix A = Matrix(2, 3);
    Matrix B = Matrix(3, 4);

    for (int i = 0; i < A.row; i++){
        for (int j = 0; j < A.col; j++){
            A[i][j] = rand() % 6 - 3;
        }
    }

    for (int i = 0; i < B.row; i++){
        for (int j = 0; j < B.col; j++){
            B[i][j] = rand() % 6 - 3;
        }
    }

    Matrix C = A * B;
    assert(verify(A, B, C));

    for (int i = 0; i < C.row; i++){
        for (int j = 0; j < C.col; j++){
            C[i][j]--;
            assert(!verify(A, B, C));
            C[i][j] += 2;
            assert(!verify(A, B, C));
            C[i][j]--;
        }
    }

    Matrix P(2, 2), Q(2, 2), R(2, 2);
    P[0][0] = 1, P[0][1] = 2, P[1][0] = 3, P[1][1] = 4;
    Q[0][0] = 5, Q[0][1] = 6, Q[1][0] = 7, Q[1][1] = 8;
    R[0][0] = 19, R[0][1] = 22, R[1][0] = 43, R[1][1] = 50;
    assert(P * Q == R && verify(P, Q, R));
    assert(Matrix(2, 3, 1) * Matrix(3, 3, 7) == Matrix(2, 3, 7));

    return 0;
}
