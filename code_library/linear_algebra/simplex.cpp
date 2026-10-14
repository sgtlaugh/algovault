/***
 *
 * Simplex
 * Maximizes or minimizes a linear objective over linear constraints on variables x_j >= 0
 *
 * Complexity: O(n m) time per pivot and O(n m) memory, the pivot count is exponential in the worst case
 *
 * Simplex lp(objective, sense): sense is Simplex::MAXIMIZE or Simplex::MINIMIZE, objective[j] is the coefficient of x_j
 * lp.add_constraint(coef, lim, cmp): coef . x <= lim, >= lim or == lim for cmp Simplex::LESSEQ, GREATEQ or EQUAL, coef has n entries
 * lp.solve(): Simplex::FEASIBLE, INFEASIBLE or UNBOUNDED, asserts it is called once per instance
 * lp.value: the optimum after FEASIBLE
 * lp.val[j]: the value of x_j at the optimum
 *
 * A free variable x becomes x1 - x2 with x1, x2 >= 0, and |x| <= M becomes the two rows x <= M and -x <= M
 * A ratio (c . x + c0) / (d . x + d0) with a positive denominator becomes linear with t = 1 / (d . x + d0) and y = t x:
 * optimize c . y + c0 t subject to d . y + d0 t == 1, every row a . x <= b turned into a . y - b t <= 0, then x = y / t
 * Each row is divided by its largest coefficient, so rows of very different size mix safely under the absolute EPS = 1e-9
 * An LP whose answer flips when a limit moves by about 1e-7 of its row, like nearly parallel rows, is beyond double precision
 * Changing `using Float = double;` to long double sets EPS to 1e-13
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Simplex{
    using Float = double;
    static constexpr Float EPS = is_same<Float, long double>::value ? 1e-13L : 1e-9;

    /// Enums rather than ints, so passing LESSEQ as the sense or a bare number fails to compile
    enum Cmp{LESSEQ, EQUAL, GREATEQ};
    enum Sense{MINIMIZE = -1, MAXIMIZE = 1};
    enum Status{INFEASIBLE, FEASIBLE, UNBOUNDED};

    Float value = 0;
    vector<Float> val;
    int n, m;
    Sense sense;
    vector<vector<Float>> ar;
    vector<Float> rhs;
    vector<int> nonzero, nonbasic, basic;

    Simplex(const vector<Float>& objective, Sense sense)
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
        /// Dantzig's rule, falling back to Bland's after a degenerate pivot since only degenerate pivots can cycle
        for (bool bland = false; ; ){
            int x = 0, y = -1;
            Float v = 0, mx = EPS;
            for (int j = 0; j <= n; j++){
                if (ar[z][j] > mx && (!bland || y < 0 || nonbasic[j] < nonbasic[y])) mx = bland ? EPS : ar[z][j], y = j;
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
            if (abs(row[j]) > EPS) nonzero[len++] = j;
        }

        /// Row 0 is the objective and row m + 1 the phase 1 objective while it exists, their rhs holds minus their value
        for (int i = 0; i < (int)ar.size(); i++){
            vector<Float>& cur = ar[i];
            if (i == x || abs(cur[y]) <= EPS) continue;
            rhs[i] -= cur[y] * rhs[x], v = cur[y], cur[y] = 0;
            for (int k = 0; k < len; k++) cur[nonzero[k]] -= v * row[nonzero[k]];
        }
    }
};

int main(){
    /// Maximize 3x0 + 2x1 subject to x0 + x1 <= 4 and 2x0 + x1 <= 5, best at the corner (1, 3)
    Simplex lp({3, 2}, Simplex::MAXIMIZE);
    lp.add_constraint({1, 1}, 4, Simplex::LESSEQ);
    lp.add_constraint({2, 1}, 5, Simplex::LESSEQ);
    assert(lp.solve() == Simplex::FEASIBLE);
    assert(abs(lp.value - 9) < 1e-9);
    assert(abs(lp.val[0] - 1) < 1e-9);
    assert(abs(lp.val[1] - 3) < 1e-9);

    /// Minimize x0 + 2x1 subject to x0 + x1 >= 3 and x0 - x1 == 1, so x0 = x1 + 1 and the cost 3x1 + 1 is least at x1 = 1
    Simplex diet({1, 2}, Simplex::MINIMIZE);
    diet.add_constraint({1, 1}, 3, Simplex::GREATEQ);
    diet.add_constraint({1, -1}, 1, Simplex::EQUAL);
    assert(diet.solve() == Simplex::FEASIBLE);
    assert(abs(diet.value - 4) < 1e-9);
    assert(abs(diet.val[0] - 2) < 1e-9);
    assert(abs(diet.val[1] - 1) < 1e-9);

    Simplex contradiction({1, 1}, Simplex::MAXIMIZE);
    contradiction.add_constraint({1, 1}, 10, Simplex::LESSEQ);
    contradiction.add_constraint({1, 1}, 20, Simplex::GREATEQ);
    assert(contradiction.solve() == Simplex::INFEASIBLE);

    Simplex runaway({1, 0}, Simplex::MAXIMIZE);
    runaway.add_constraint({1, -1}, 1, Simplex::LESSEQ);
    assert(runaway.solve() == Simplex::UNBOUNDED);

    return 0;
}
