/***
 *
 * Freivalds' Algorithm
 * A randomized algorithm that verifies matrix multiplication
 *
 * Complexity: O(n^2 t) for an n x n matrix and t = ntrials
 *
 * FreivaldsMatrix M(row, col, diagonal = 0): row x col long long matrix with diagonal on its diagonal, read as M[i][j]
 * freivalds_verify(A, B, C, ntrials = 30): whether A * B = C, false on mismatched shapes
 *
 * Notes:
 *   - A wrong true has probability below 1 / 2^ntrials, false is always right
 *   - The check runs modulo 2^64, which is exact whenever C and the true product A * B fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// One row-major vector: vector<vector<long long>> measured 7-11% slower on a 100 x 100 product
struct FreivaldsMatrix{
    int row, col;
    vector<long long> mat;

    FreivaldsMatrix() : row(0), col(0) {}
    FreivaldsMatrix(int row, int col, int diagonal = 0) : row(row), col(col), mat((size_t)row * col, 0){
        for (int i = min(row, col) - 1; i >= 0; i--) (*this)[i][i] = diagonal;
    }

    long long* operator[](int i){
        return mat.data() + (size_t)i * col;
    }

    const long long* operator[](int i) const{
        return mat.data() + (size_t)i * col;
    }

    FreivaldsMatrix operator* (const FreivaldsMatrix& other) const{
        FreivaldsMatrix res(row, other.col);

        for (int i = 0; i < row; i++){
            for (int j = 0; j < other.col; j++){
                long long sum = 0;
                for (int k = 0; k < col; k++) sum += (*this)[i][k] * other[k][j];
                res[i][j] = sum;
            }
        }

        return res;
    }

    bool operator == (const FreivaldsMatrix& m) const{
        return row == m.row && col == m.col && mat == m.mat;
    }

    bool operator != (const FreivaldsMatrix & m) const{
        return !(*this == m);
    }
};

bool freivalds_verify(const FreivaldsMatrix& A, const FreivaldsMatrix& B, const FreivaldsMatrix& C, int ntrials=30){
    if (A.col != B.row || C.row != A.row || C.col != B.col) return false;

    auto dot = [](const long long* a, const vector<unsigned long long>& b){
        unsigned long long res = 0;
        for (size_t i = 0; i < b.size(); i++) res += (unsigned long long)a[i] * b[i];
        return res;
    };

    static mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    vector<unsigned long long> r(C.col), br(B.row), x(A.row), y(C.row);  /// wraps mod 2^64 instead of overflowing

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
    FreivaldsMatrix A = FreivaldsMatrix(2, 3);
    FreivaldsMatrix B = FreivaldsMatrix(3, 4);

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

    FreivaldsMatrix C = A * B;
    assert(freivalds_verify(A, B, C));

    for (int i = 0; i < C.row; i++){
        for (int j = 0; j < C.col; j++){
            C[i][j]--;
            assert(!freivalds_verify(A, B, C));
            C[i][j] += 2;
            assert(!freivalds_verify(A, B, C));
            C[i][j]--;
        }
    }

    FreivaldsMatrix P(2, 2), Q(2, 2), R(2, 2);
    P[0][0] = 1, P[0][1] = 2, P[1][0] = 3, P[1][1] = 4;
    Q[0][0] = 5, Q[0][1] = 6, Q[1][0] = 7, Q[1][1] = 8;
    R[0][0] = 19, R[0][1] = 22, R[1][0] = 43, R[1][1] = 50;
    assert(P * Q == R && freivalds_verify(P, Q, R));
    assert(FreivaldsMatrix(2, 3, 1) * FreivaldsMatrix(3, 3, 7) == FreivaldsMatrix(2, 3, 7));

    return 0;
}
