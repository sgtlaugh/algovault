/***
 *
 * Thomas Algorithm
 * Solves a tridiagonal system of n linear equations in linear time
 *
 * Complexity: O(n) time and memory
 *
 * solve_tridiagonal(lower, diagonal, upper, rhs): x with lower[i] x[i-1] + diagonal[i] x[i] + upper[i] x[i+1] = rhs[i]
 * solve_tridiagonal_m_matrix(lower, upper, excess, rhs): x for diagonals -lower, lower + upper + excess, -upper
 *
 * Notes:
 *   - All vectors have size n and a floating type, lower[0] and upper[n - 1] lie outside the matrix and multiply no x
 *   - Pivots stay nonzero when every |diagonal[i]| > |lower[i]| + |upper[i]|, taking lower[0] and upper[n - 1] as 0
 *   - Weak dominance (>=) does too when row 0 is strict and every lower[i] past the first is nonzero, like -1 2 -1
 *   - Otherwise a pivot can reach 0 and the answer turns inf or nan
 *   - Nonzero pivots do not mean accuracy: hitting probabilities of a walk with strong drift are too ill-conditioned
 *   - solve_tridiagonal_m_matrix needs every input >= 0 and never subtracts, so every x[i] keeps a small relative error
 *   - In it lower[0] and upper[n - 1] count in the diagonal as steps out of the matrix, excess[i] is any further leak
 *   - Take excess from the model, like 0 or a leak chance: computing diagonal - lower - upper brings cancellation back
 *   - A block of rows with no excess and no step out of the matrix makes it singular, and the answer turns inf or nan
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

template <class T>
vector<T> solve_tridiagonal_m_matrix(const vector<T>& lower, const vector<T>& upper, const vector<T>& excess, vector<T> rhs){
    int n = rhs.size();
    assert((int)lower.size() == n && (int)upper.size() == n && (int)excess.size() == n);
    vector<T> ratio(n);

    /// Forward sweep with pivot upper[i] + surplus, surplus = diagonal - lower[i] * ratio[i - 1] by sums and products
    T surplus = 0, inverse_pivot = 1;
    for (int i = 0; i < n; i++){
        /// lower[i] (1 - ratio[i - 1]) without subtracting; all of lower[0] leaves the matrix, row 0 has no row above
        T carried = i ? lower[i] * surplus * inverse_pivot : lower[i];
        surplus = excess[i] + carried;
        inverse_pivot = 1 / (upper[i] + surplus);
        ratio[i] = upper[i] * inverse_pivot;
        rhs[i] = (rhs[i] + (i ? lower[i] * rhs[i - 1] : 0)) * inverse_pivot;
    }

    /// Back substitution, rhs turns into x from the last row up
    for (int i = n - 2; i >= 0; i--) rhs[i] += ratio[i] * rhs[i + 1];

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

    /// Fair walk on 0..4 stopping at both ends reaches 4 first from i with chance i / 4, the step 3 -> 4 is the rhs 1/2
    vector<double> half = {0.5, 0.5, 0.5}, none = {0, 0, 0}, reach = {0, 0, 0.5};
    vector<double> h = solve_tridiagonal_m_matrix(half, half, none, reach);
    assert(h.size() == 3);
    assert(abs(h[0] - 0.25) < 1e-12);
    assert(abs(h[1] - 0.5) < 1e-12);
    assert(abs(h[2] - 0.75) < 1e-12);

    return 0;
}
