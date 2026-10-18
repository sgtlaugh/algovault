/***
 *
 * Thomas Algorithm
 * Solves a tridiagonal system of n linear equations in linear time
 *
 * Complexity: O(n) time and memory
 *
 * solve_tridiagonal(lower, diagonal, upper, rhs): x with lower[i] x[i - 1] + diagonal[i] x[i] + upper[i] x[i + 1] = rhs[i] for every i
 * solve_tridiagonal_m_matrix(lower, upper, excess, rhs): x with (lower[i] + upper[i] + excess[i]) x[i] - lower[i] x[i - 1] - upper[i] x[i + 1] = rhs[i], all inputs >= 0
 *
 * All four vectors have size n and a floating type, lower[0] and upper[n - 1] lie outside the matrix and multiply no x
 * No pivoting: every pivot stays nonzero when every |diagonal[i]| > |lower[i]| + |upper[i]|, entries outside the matrix counting as 0
 * Weak dominance (>=) keeps the pivots nonzero too when the first row is strict and every lower[i] past the first is nonzero, like -1 2 -1
 * Otherwise a pivot can reach 0 and the answer turns inf or nan
 * Nonzero pivots do not mean an accurate answer: hitting probabilities of a random walk with strong drift are weakly dominant yet too ill-conditioned for floating point
 * solve_tridiagonal_m_matrix never subtracts, so every x[i] keeps a small relative error however ill-conditioned the system, as for those random walks
 * Its lower[0] and upper[n - 1] still count in the diagonal, as steps out of the matrix, and excess[i] is any further leak from row i
 * Take excess from the model, like 0 or a leak probability: computing it as diagonal - lower - upper brings the cancellation back
 * A block of rows with no excess and no step out of the matrix makes it singular, and the answer turns inf or nan
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

    /// Forward sweep as in solve_tridiagonal, but the pivot is upper[i] plus a surplus built only from sums and products
    /// The surplus is diagonal - lower[i] * ratio[i - 1], and lower[i] (1 - ratio[i - 1]) = lower[i] * surplus[i - 1] / pivot[i - 1]
    T surplus = 0, inverse_pivot = 1;
    for (int i = 0; i < n; i++){
        T carried = i ? lower[i] * surplus * inverse_pivot : lower[i];  /// row 0 has no row above, all of lower[0] leaves the matrix
        surplus = excess[i] + carried;
        inverse_pivot = 1 / (upper[i] + surplus);
        ratio[i] = upper[i] * inverse_pivot;
        rhs[i] = (rhs[i] + (i ? lower[i] * rhs[i - 1] : 0)) * inverse_pivot;
    }

    /// Back substitution, rhs turns into x from the last row up
    for (int i = n - 2; i >= 0; i--) rhs[i] += ratio[i] * rhs[i + 1];

    return rhs;
}

/// Chance to reach hi before lo from each cell lo + 1..hi - 1 of the line, p[i] = chance to pass from i to i + 1
vector<double> reach_high_first(const vector<double>& p, int lo, int hi){
    vector<double> lower, upper, excess, rhs;
    for (int i = lo + 1; i < hi; i++){
        lower.push_back(1 - p[i]);
        upper.push_back(p[i]);
        excess.push_back(0);
        rhs.push_back(0);
    }
    rhs.back() = upper.back();  /// the step into hi carries the value 1
    return solve_tridiagonal_m_matrix(lower, upper, excess, rhs);
}

int main(){
    int n, k;
    while (scanf("%d %d", &n, &k) == 2){
        vector<double> p(n, 0);
        for (int i = 1; i < n; i++) scanf("%lf", &p[i]);

        /// Person n wins when the token covers 1 and n - 1 before stepping onto n, cutting the circle at n leaves a line 0..n
        if (n == 2){
            puts("1.0000000000");
            continue;
        }

        /// first_left: from k, the token meets 1 before n - 1
        double first_left = k == 1 ? 1 : k == n - 1 ? 0 : 1 - reach_high_first(p, 1, n - 1)[k - 2];

        /// From 1 it must then reach n - 1 before stepping left onto n, which sits at 0 on the line
        double left_then_right = reach_high_first(p, 0, n - 1)[0];

        /// From n - 1 it must reach 1 before stepping right onto n, the same walk seen in a mirror
        vector<double> mirrored(n, 0);
        for (int i = 1; i < n; i++) mirrored[n - i] = 1 - p[i];
        double right_then_left = reach_high_first(mirrored, 0, n - 1)[0];

        printf("%.10f\n", first_left * left_then_right + (1 - first_left) * right_then_left);
    }
    return 0;
}
