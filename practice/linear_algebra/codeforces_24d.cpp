/***
 *
 * Thomas Algorithm
 * Solves a tridiagonal system of n linear equations in linear time
 *
 * Complexity: O(n) time and memory
 *
 * solve_tridiagonal(lower, diagonal, upper, rhs): x with lower[i] x[i - 1] + diagonal[i] x[i] + upper[i] x[i + 1] = rhs[i] for every i
 *
 * All four vectors have size n and a floating type, lower[0] and upper[n - 1] lie outside the matrix and are ignored
 * No pivoting: every pivot stays nonzero when every |diagonal[i]| > |lower[i]| + |upper[i]|, entries outside the matrix counting as 0
 * Weak dominance (>=) keeps the pivots nonzero too when the first row is strict and every lower[i] past the first is nonzero, like -1 2 -1
 * Otherwise a pivot can reach 0 and the answer turns inf or nan
 * Nonzero pivots do not mean an accurate answer: hitting probabilities of a random walk with strong drift are weakly dominant yet too ill-conditioned for floating point
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
vector<T> solve_tridiagonal(vector<T> lower, vector<T> diagonal, vector<T> upper, vector<T> rhs){
    int n = diagonal.size();
    assert((int)lower.size() == n && (int)upper.size() == n && (int)rhs.size() == n);
    if (n == 0) return {};

    /// Forward sweep, row i becomes x[i] + upper[i] x[i + 1] = rhs[i]
    upper[0] /= diagonal[0];
    rhs[0] /= diagonal[0];
    for (int i = 1; i < n; i++){
        T scale = 1 / (diagonal[i] - lower[i] * upper[i - 1]);
        upper[i] *= scale;
        rhs[i] = (rhs[i] - lower[i] * rhs[i - 1]) * scale;
    }

    /// Back substitution, rhs turns into x from the last row up
    for (int i = n - 2; i >= 0; i--) rhs[i] -= upper[i] * rhs[i + 1];

    return rhs;
}

int main(){
    int rows, cols, start_row, start_col;
    scanf("%d %d %d %d", &rows, &cols, &start_row, &start_col);

    /// The robot never moves up, so each row only needs the expectations of the row below, solved bottom up
    vector<double> below(cols, 0);
    for (int r = rows - 1; r >= start_row; r--){
        vector<double> lower(cols), diagonal(cols), upper(cols), rhs(cols);
        for (int c = 0; c < cols; c++){
            int moves = 2 + (c > 0) + (c + 1 < cols);
            diagonal[c] = 1 - 1.0 / moves;
            lower[c] = c > 0 ? -1.0 / moves : 0;
            upper[c] = c + 1 < cols ? -1.0 / moves : 0;
            rhs[c] = 1 + below[c] / moves;
        }
        below = solve_tridiagonal(lower, diagonal, upper, rhs);
    }

    printf("%.10f\n", below[start_col - 1]);
    return 0;
}
