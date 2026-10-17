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
    /// 2x + y = 4, x + 3y + z = 9, y + 2z = 6 has the solution x = 1, y = 2, z = 2
    vector<double> lower = {0, 1, 1}, diagonal = {2, 3, 2}, upper = {1, 1, 0}, rhs = {4, 9, 6};
    vector<double> x = solve_tridiagonal(lower, diagonal, upper, rhs);
    assert(x.size() == 3);
    assert(abs(x[0] - 1) < 1e-12);
    assert(abs(x[1] - 2) < 1e-12);
    assert(abs(x[2] - 2) < 1e-12);

    return 0;
}
