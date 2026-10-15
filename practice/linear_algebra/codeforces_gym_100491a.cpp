/***
 *
 * Simplex
 * Maximizes or minimizes a linear objective subject to linear constraints on variables x_j >= 0
 *
 * Complexity: O(n m) time per pivot and O(n m) memory, the pivot count is exponential in the worst case
 *
 * Simplex lp(objective, sense): an LP in n = objective.size() variables, sense is Simplex::MAXIMIZE or Simplex::MINIMIZE
 * BasicSimplex<long double> lp(objective, sense): the same in long double, about 4x slower, resolving 1e-13 instead of 1e-9
 * lp.add_constraint(coef, lim, cmp): adds coef . x <= lim, >= lim or == lim for cmp Simplex::LESSEQ, GREATEQ or EQUAL, coef has n entries
 * lp.solve(): Simplex::FEASIBLE, INFEASIBLE or UNBOUNDED, asserts it is called once per instance
 * lp.value: the optimum after FEASIBLE
 * lp.val[j]: the value of x_j at the optimum
 *
 * A free variable x becomes x1 - x2 with x1, x2 >= 0
 * |x| <= M becomes the two rows x <= M and -x <= M
 * A ratio objective with a denominator positive on every feasible x becomes linear with t = 1 / denominator and y = t x:
 *
 *     max (3x1 + 2x2 + 4x3 + 6) / (3x1 + 3x2 + 2x3 + 5)    ->    max 3y1 + 2y2 + 4y3 + 6t
 *     2x1 + 3x2 + 5x3 >= 23                                      3y1 + 3y2 + 2y3 + 5t == 1
 *     3x1 + 5x2 + 4x3 <= 30                                      2y1 + 3y2 + 5y3 - 23t >= 0
 *                                                                3y1 + 5y2 + 4y3 - 30t <= 0
 *
 *   then x = y / t
 *
 * Rows are scaled internally, so mixing very large and very small coefficients is fine
 * The answer is unreliable only when moving a limit by about 1e-7 of its row changes it, as with nearly parallel rows
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Enums rather than ints, so passing LESSEQ as the sense or a bare number fails to compile
/// They live outside the template so Simplex::MAXIMIZE also serves BasicSimplex<long double>
struct SimplexOptions{
    enum Cmp{LESSEQ, EQUAL, GREATEQ};
    enum Sense{MINIMIZE = -1, MAXIMIZE = 1};
    enum Status{INFEASIBLE, FEASIBLE, UNBOUNDED};
};

template <typename Float>
struct BasicSimplex : SimplexOptions{
    static constexpr bool WIDE = is_same<Float, long double>::value;
    /// EPS decides pivots, feasibility and optimality; elimination drops only entries below DROP, since dropping at EPS
    /// loses entries that still matter once variables grow large
    static constexpr Float EPS = WIDE ? 1e-13L : 1e-9, DROP = WIDE ? 1e-16L : 1e-12;

    Float value = 0;
    vector<Float> val;
    int n, m;
    Sense sense;
    vector<vector<Float>> ar;
    vector<Float> rhs, norm;
    vector<int> nonzero, nonbasic, basic;

    BasicSimplex(const vector<Float>& objective, Sense sense)
        : val(objective.size(), 0), n(objective.size()), m(0), sense(sense), ar(1, vector<Float>(n + 1, 0)), rhs(1, 0), nonzero(n + 1), nonbasic(n + 1){
        for (int j = 0; j < n; j++) ar[0][j + 1] = objective[j] * sense;
    }

    void add_constraint(const vector<Float>& coef, Float lim, Cmp cmp){
        assert((int)coef.size() == n);
        if (cmp != GREATEQ) add_row(coef, lim, 1);
        if (cmp != LESSEQ) add_row(coef, lim, -1);
    }

    Status solve(){
        assert(basic.empty());  /// the tableau is pivoted in place, a second call would answer from garbage

        /// Variable ids: 0 is the phase 1 artificial, 1..n the inputs, n + i the slack of row i
        basic.assign(m + 1, 0);
        for (int j = 0; j <= n; j++) nonbasic[j] = j;
        for (int i = 1; i <= m; i++) basic[i] = i + n;

        /// Phase 1: relax every row by an artificial x0 in column 0, then maximize -x0 to reach a feasible basis
        int x = 1;
        for (int i = 2; i <= m; i++){
            if (rhs[i] < rhs[x]) x = i;
        }
        if (m > 0 && rhs[x] < -EPS){
            ar.emplace_back(n + 1, 0), rhs.push_back(0);
            for (int i = 1; i <= m + 1; i++) ar[i][0] = -1;

            pivot(x, 0);
            optimize(m + 1);
            if (rhs[m + 1] > EPS) return INFEASIBLE;

            /// x0 may stay basic at a value within EPS and grow in phase 2, so pivot it out, a row with no pivot is redundant
            for (int i = 1; i <= m; i++){
                if (basic[i] != 0) continue;
                int y = 0;
                for (int j = 1; j <= n; j++){
                    if (abs(ar[i][j]) > abs(ar[i][y])) y = j;
                }
                if (abs(ar[i][y]) > EPS) pivot(i, y);
            }
            for (int j = 0; j <= n; j++){
                if (nonbasic[j] == 0) for (auto& row : ar) row[j] = 0;
            }
            ar.pop_back(), rhs.pop_back();
        }

        if (!optimize(0)) return UNBOUNDED;

        value = -rhs[0] * sense;
        vector<int> row_of(n + 1, 0);
        for (int i = 1; i <= m; i++){
            if (basic[i] <= n) row_of[basic[i]] = i;
        }
        for (int j = 1; j <= n; j++) val[j - 1] = row_of[j] ? rhs[row_of[j]] : 0;
        return FEASIBLE;
    }

    void add_row(const vector<Float>& coef, Float lim, int sign){
        /// Dividing by the largest coefficient makes the absolute EPS mean the same in every row
        Float scale = 0;
        for (Float c : coef) scale = max(scale, abs(c));
        if (scale == 0) scale = 1;
        m++;
        rhs.push_back(lim * sign / scale);
        ar.emplace_back(n + 1, 0);
        for (int j = 0; j < n; j++) ar[m][j + 1] = coef[j] * sign / scale;
    }

    bool optimize(int z){
        /// Steepest edge, falling back to Bland's rule after a degenerate pivot since only degenerate pivots can cycle
        for (bool bland = false; ; ){
            int x = 0, y = -1;
            Float v = 0;
            if (bland){
                for (int j = 0; j <= n; j++){
                    if (ar[z][j] > EPS && (y < 0 || nonbasic[j] < nonbasic[y])) y = j;
                }
            }
            else{
                /// Gain per unit length of the edge, not per unit of the variable: 1.4-3.6x faster from n = 300, 20-30% slower under 15 ms
                norm.assign(n + 1, 1);
                for (int i = 1; i <= m; i++){
                    for (int j = 0; j <= n; j++) norm[j] += ar[i][j] * ar[i][j];
                }
                Float best = 0;
                for (int j = 0; j <= n; j++){
                    if (ar[z][j] > EPS && ar[z][j] * ar[z][j] / norm[j] > best) best = ar[z][j] * ar[z][j] / norm[j], y = j;
                }
            }
            if (y < 0) return true;

            for (int i = 1; i <= m; i++){
                if (ar[i][y] <= EPS) continue;
                Float u = rhs[i] / ar[i][y];
                /// Near ties take the largest pivot to limit rounding, the lowest id only while Bland's rule prevents cycling
                bool better = bland ? basic[i] < basic[x] : ar[i][y] > ar[x][y];
                if (x == 0 || u < v - EPS || (u <= v + EPS && better)) x = i, v = u;
            }
            if (x == 0) return false;
            bland = v <= EPS;
            pivot(x, y);
        }
    }

    void pivot(int x, int y){
        int len = 0;
        vector<Float>& row = ar[x];
        Float v = row[y];

        swap(basic[x], nonbasic[y]);
        rhs[x] /= v, row[y] = 1;
        for (int j = 0; j <= n; j++){
            row[j] /= v;
            if (abs(row[j]) > DROP) nonzero[len++] = j;
        }

        /// Row 0 is the objective and row m + 1 the phase 1 objective while it exists, their rhs holds minus their value
        for (int i = 0; i < (int)ar.size(); i++){
            vector<Float>& cur = ar[i];
            if (i == x || abs(cur[y]) <= DROP) continue;
            rhs[i] -= cur[y] * rhs[x], v = cur[y], cur[y] = 0;
            for (int k = 0; k < len; k++) cur[nonzero[k]] -= v * row[nonzero[k]];
        }
    }
};

using Simplex = BasicSimplex<double>;

int main(){
    freopen("arbitrage.in", "r", stdin);
    freopen("arbitrage.out", "w", stdout);

    int n, m;
    scanf("%d %d", &n, &m);
    vector<int> from(m), to(m), cap(m);
    vector<double> profit(m);
    for (int e = 0; e < m; e++){
        int a;
        scanf("%d %d %d %d", &from[e], &to[e], &a, &cap[e]);
        profit[e] = a;
    }

    /// A circulation: every currency nets to zero; the incidence matrix is totally unimodular, so the optimum is integral
    Simplex lp(profit, Simplex::MAXIMIZE);
    for (int v = 1; v <= n; v++){
        vector<double> net(m, 0);
        for (int e = 0; e < m; e++){
            if (from[e] == v) net[e] += 1;
            if (to[e] == v) net[e] -= 1;
        }
        lp.add_constraint(net, 0, Simplex::EQUAL);
    }
    for (int e = 0; e < m; e++){
        vector<double> one(m, 0);
        one[e] = 1;
        lp.add_constraint(one, cap[e], Simplex::LESSEQ);
    }

    lp.solve();
    printf("%lld\n", llround(lp.value));
    return 0;
}
