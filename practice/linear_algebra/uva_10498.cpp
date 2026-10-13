/***
 *
 * Simplex Algorithm for Linear Programming
 * Solves maximize/minimize linear objective subject to linear constraints
 *
 * Complexity: O(n * m) per pivot and O(n * m) memory, the pivot count is exponential in the worst case
 * Optimized with sparse pivot - only updates non-zero coefficients
 *
 * Usage:
 *   1. Simplex lp(objective, Simplex::MAXIMIZE or Simplex::MINIMIZE), objective[j] is the coefficient of x_j
 *   2. lp.add_constraint(coefficients, rhs, Simplex::LESSEQ or Simplex::GREATEQ or Simplex::EQUAL)
 *   3. lp.solve() returns Simplex::FEASIBLE, Simplex::INFEASIBLE or Simplex::UNBOUNDED
 *   4. If feasible, lp.value is the optimum and lp.val[j] the value of x_j
 *
 * Variables x_0, ..., x_{n - 1} satisfy x_j >= 0, n = objective.size() and every coefficient vector has n entries
 * solve() pivots the tableau in place, so it asserts it is called once per instance
 * For unbounded variable: replace x with (x1 - x2) where x1, x2 >= 0
 * For constraint |x| <= M: add both x <= M and -x <= M
 *
 * Change `using Float = double;` to long double for higher precision (EPS follows)
 *
 * Example:
 *   Maximize 3x0 + 2x1 subject to x0 + x1 <= 4, 2x0 + x1 <= 5
 *
 *   Simplex lp({3, 2}, Simplex::MAXIMIZE);
 *   lp.add_constraint({1, 1}, 4, Simplex::LESSEQ);
 *   lp.add_constraint({2, 1}, 5, Simplex::LESSEQ);
 *   if (lp.solve() == Simplex::FEASIBLE) ...  // lp.value = 9, lp.val = {1, 3}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Simplex{
    using Float = double;
    static constexpr Float EPS = is_same<Float, long double>::value ? 1e-13L : 1e-9;

    static constexpr int LESSEQ = -1, EQUAL = 0, GREATEQ = 1;
    static constexpr int MINIMIZE = -1, MAXIMIZE = 1;
    static constexpr int INFEASIBLE = -1, FEASIBLE = 1, UNBOUNDED = 666;

    Float value = 0;
    vector<Float> val;

    Simplex(const vector<Float>& objective, int sense)
        : val(objective.size(), 0), n(objective.size()), m(0), flag(sense), ar(1, vector<Float>(n + 1, 0)), rhs(1, 0), adj(n + 1), down(n + 1){
        for (int j = 0; j < n; j++) ar[0][j + 1] = objective[j] * flag;
    }

    void add_constraint(const vector<Float>& coef, Float lim, int cmp){
        assert((int)coef.size() == n);
        if (cmp != GREATEQ) add_row(coef, lim, 1);
        if (cmp != LESSEQ) add_row(coef, lim, -1);
    }

    int solve(){
        /// link is filled only by solve, a second call would read the already pivoted tableau and return garbage
        assert(link.empty());
        int i, j, x = 1;

        /// Variable ids: 0 is the phase 1 artificial, 1..n the inputs, n + i the slack of row i
        link.assign(m + 1, 0);
        for (i = 0; i <= n; i++) down[i] = i;
        for (i = 1; i <= m; i++) link[i] = i + n;

        // Phase 1: relax every row by an artificial x0 in column 0, then maximize -x0 to reach a feasible basis
        for (i = 2; i <= m; i++){
            if (rhs[i] < rhs[x]) x = i;
        }
        if (m > 0 && rhs[x] < -EPS){
            ar.emplace_back(n + 1, 0), rhs.push_back(0);
            for (i = 1; i <= m; i++) ar[i][0] = -1;
            ar[m + 1][0] = -1;

            pivot(x, 0);
            optimize(m + 1);
            if (rhs[m + 1] > EPS) return INFEASIBLE;

            /// Ratio ties go to x0's id 0, but x0 can stay basic at a value within EPS and would grow in phase 2, so pivot it out
            /// A row with no usable pivot is redundant and stays inert
            for (i = 1; i <= m; i++){
                if (link[i] != 0) continue;
                int y = 0;
                for (j = 1; j <= n; j++){
                    if (abs(ar[i][j]) > abs(ar[i][y])) y = j;
                }
                if (abs(ar[i][y]) > EPS) pivot(i, y);
            }
            for (j = 0; j <= n; j++){
                if (down[j] == 0) for (auto& r : ar) r[j] = 0;
            }
            ar.pop_back(), rhs.pop_back();
        }

        // Phase 2: optimize objective function
        if (!optimize(0)) return UNBOUNDED;

        value = -rhs[0] * flag;
        vector<int> idx(n + 1, 0);
        for (i = 1; i <= m; i++){
            if (link[i] <= n) idx[link[i]] = i;
        }

        for (j = 1; j <= n; j++) val[j - 1] = idx[j] ? rhs[idx[j]] : 0;
        return FEASIBLE;
    }

    int n, m, flag;
    vector<vector<Float>> ar;
    vector<Float> rhs;
    vector<int> adj, down, link;

    void add_row(const vector<Float>& coef, Float lim, int sign){
        m++;
        rhs.push_back(lim * sign);
        ar.emplace_back(n + 1, 0);
        for (int j = 0; j < n; j++) ar[m][j + 1] = coef[j] * sign;
    }

    bool optimize(int z){
        /// Dantzig's rule, falling back to Bland's after a degenerate pivot since only degenerate pivots can cycle
        for (bool bland = false; ; ){
            int i, x = 0, y = -1;
            Float u, v = 0, mx = EPS;
            for (i = 0; i <= n; i++){
                if (ar[z][i] > mx && (!bland || y < 0 || down[i] < down[y])) mx = bland ? EPS : ar[z][i], y = i;
            }
            if (y < 0) return true;

            for (i = 1; i <= m; i++){
                if (ar[i][y] > EPS){
                    u = rhs[i] / ar[i][y];
                    if (x == 0 || u < v - EPS || (u <= v + EPS && link[i] < link[x])) x = i, v = u;
                }
            }
            if (x == 0) return false;
            bland = v <= EPS;
            pivot(x, y);
        }
    }

    void pivot(int x, int y){
        int i, j, len = 0;
        vector<Float>& row = ar[x];
        Float v = row[y];

        swap(link[x], down[y]);
        rhs[x] /= v, row[y] = 1;
        for (j = 0; j <= n; j++){
            row[j] /= v;
            if (abs(row[j]) > EPS) adj[len++] = j;
        }

        /// Row 0 is the objective and row m + 1 the phase 1 objective while it exists, their rhs holds minus their value
        for (i = 0; i < (int)ar.size(); i++){
            vector<Float>& cur = ar[i];
            if (abs(cur[y]) > EPS && i != x){
                rhs[i] -= cur[y] * rhs[x], v = cur[y], cur[y] = 0;
                for (j = 0; j < len; j++) cur[adj[j]] -= (v * row[adj[j]]);
            }
        }
    }
};

int main(){
    int n, m;
    while (scanf("%d %d", &n, &m) == 2){
        /// Every contestant gets the same amounts x, so the bill is m times the price of one share
        vector<double> price(n), happiness(n);
        for (auto& p : price) scanf("%lf", &p), p *= m;
        Simplex lp(price, Simplex::MAXIMIZE);

        for (int i = 0; i < m; i++){
            double total;
            for (auto& h : happiness) scanf("%lf", &h);
            scanf("%lf", &total);
            lp.add_constraint(happiness, total, Simplex::LESSEQ);
        }
        lp.solve();
        printf("Nasa can spend %d taka.\n", (int)ceil(lp.value));
    }
    return 0;
}
