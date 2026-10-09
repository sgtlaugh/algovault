/***
 *
 * Simplex Algorithm for Linear Programming
 * Solves maximize/minimize linear objective subject to linear constraints
 *
 * Complexity: O(2^n) worst case, polynomial average case, O(n * m) per pivot and O(n * m) memory
 * Optimized with sparse pivot - only updates non-zero coefficients
 *
 * Usage:
 *   1. Simplex lp(n, objective, MAXIMIZE or MINIMIZE)
 *   2. lp.add_constraint(coefficients, rhs, LESSEQ or GREATEQ or EQUAL)
 *   3. lp.solve(result) returns FEASIBLE, INFEASIBLE, or UNBOUNDED
 *   4. If feasible, result has optimal value, lp.val[] has variable values
 *
 * Variables x1, x2, ..., xn (1-indexed) satisfy xi >= 0
 * objective and coefficients have n + 1 entries, index 0 is ignored
 * solve() pivots the tableau in place, so it asserts it is called once per instance
 * Phase 1 pivots on a random column without a ratio test, so it carries no termination guarantee
 * For unbounded variable: replace x with (x1 - x2) where x1, x2 >= 0
 * For constraint |x| <= M: add both x <= M and -x <= M
 *
 * Change `using Float = long double;` for higher precision (auto-adjusts EPS)
 *
 * Example:
 *   Maximize 3x1 + 2x2 subject to x1 + x2 <= 4, 2x1 + x2 <= 5
 *
 *   Float result;
 *   Simplex lp(2, {0, 3, 2}, MAXIMIZE);
 *   lp.add_constraint({0, 1, 1}, 4, LESSEQ);
 *   lp.add_constraint({0, 2, 1}, 5, LESSEQ);
 *   int status = lp.solve(result);
 *
 *   // result = 9, lp.val[1] = 1, lp.val[2] = 3
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

using Float = double;  // Change to 'long double' for higher precision

const Float EPS = std::is_same<Float, long double>::value ? 1e-13L : 1e-9;

const int EQUAL = 0;
const int LESSEQ = -1;
const int GREATEQ = 1;
const int MINIMIZE = -1;
const int MAXIMIZE = 1;

const int FEASIBLE = 1;
const int INFEASIBLE = -1;
const int UNBOUNDED = 666;

struct Simplex{
    int n, m, flag;
    vector<vector<Float>> ar;
    vector<Float> val, rhs;
    vector<int> adj, down, link;

    Simplex(int n, const vector<Float>& obj, int min_or_max)
        : n(n), m(0), flag(min_or_max), ar(1, vector<Float>(n + 1, 0)), val(n + 1, 0), rhs(1, 0), adj(n + 1), down(n + 1){
        assert((int)obj.size() == n + 1);
        for (int i = 1; i <= n; i++) ar[0][i] = obj[i] * flag;
    }

    void add_constraint(const vector<Float>& coef, Float lim, int cmp){
        assert((int)coef.size() == n + 1);
        if (cmp != GREATEQ) add_row(coef, lim, 1);
        if (cmp != LESSEQ) add_row(coef, lim, -1);
    }

    void add_row(const vector<Float>& coef, Float lim, int sign){
        m++;
        rhs.push_back(lim * sign);
        ar.emplace_back(n + 1, 0);
        for (int i = 1; i <= n; i++) ar[m][i] = coef[i] * sign;
    }

    void pivot(int x, int y, Float& res){
        int i, j, len = 0;
        vector<Float>& row = ar[x];
        Float v = row[y];

        swap(link[x], down[y]);
        rhs[x] /= v, row[y] = 1;
        for (j = 1; j <= n; j++){
            row[j] /= v;
            if (abs(row[j]) > EPS) adj[len++] = j;
        }

        for (i = 1; i <= m; i++){
            vector<Float>& cur = ar[i];
            if (abs(cur[y]) > EPS && i != x){
                rhs[i] -= cur[y] * rhs[x], v = cur[y], cur[y] = 0;
                for (j = 0; j < len; j++) cur[adj[j]] -= (v * row[adj[j]]);
            }
        }

        vector<Float>& obj = ar[0];
        res += (obj[y] * rhs[x]), v = obj[y], obj[y] = 0;
        for (j = 0; j < len; j++) obj[adj[j]] -= (v * row[adj[j]]);
    }

    int solve(Float& res){
        /// link is filled only by solve, a second call would read the already pivoted tableau and return garbage
        assert(link.empty());
        res = 0;
        int i, x, y;
        Float u, v, mn, mx;

        link.assign(m + 1, 0);
        for (i = 1; i <= n; i++) down[i] = i;
        for (i = 1; i <= m; i++) link[i] = i + n;

        // Phase 1: find feasible solution
        while (true){
            x = 0, y = 0, mn = -EPS;
            for (i = 1; i <= m; i++){
                if (rhs[i] < mn) mn = rhs[i], x = i;
            }
            if (x == 0) break;

            for (i = 1; i <= n; i++){
                if (ar[x][i] < -EPS){
                    y = i;
                    if (rand() & 1) break;
                }
            }
            if (y == 0) return INFEASIBLE;
            pivot(x, y, res);
        }

        // Phase 2: optimize objective function
        /// Dantzig's rule, falling back to Bland's after a degenerate pivot since only degenerate pivots can cycle
        for (bool bland = false; ; ){
            x = 0, y = 0, mx = EPS;
            for (i = 1; i <= n; i++){
                if (ar[0][i] > mx && (!bland || !y || down[i] < down[y])) mx = bland ? EPS : ar[0][i], y = i;
            }
            if (y == 0) break;

            for (i = 1; i <= m; i++){
                if (ar[i][y] > EPS){
                    u = rhs[i] / ar[i][y];
                    if (x == 0 || u < v - EPS || (u <= v + EPS && link[i] < link[x])) x = i, v = u;
                }
            }
            if (x == 0) return UNBOUNDED;
            bland = v <= EPS;
            pivot(x, y, res);
        }

        res *= flag;
        vector<int> idx(n + 1, 0);
        for (i = 1; i <= m; i++){
            if (link[i] <= n) idx[link[i]] = i;
        }

        for (i = 1; i <= n; i++) val[i] = rhs[idx[i]];
        return FEASIBLE;
    }
};

int main(){
    // Production planning problem
    // Maximize profit: 5x1 + 4x2 + 3x3
    // Subject to:      2x1 + 3x2 + x3 <= 100  (material A)
    //                  4x1 + x2 + 2x3 <= 80   (material B)
    //                  3x1 + 4x2 + 2x3 <= 120 (labor hours)

    Float result = 0;

    Simplex production(3, {0, 5, 4, 3}, MAXIMIZE);
    production.add_constraint({0, 2, 3, 1}, 100, LESSEQ);
    production.add_constraint({0, 4, 1, 2}, 80, LESSEQ);
    production.add_constraint({0, 3, 4, 2}, 120, LESSEQ);

    assert(production.solve(result) == FEASIBLE);
    assert(abs(result - 153.333333333) < 1e-6);

    // Infeasible constraints
    Simplex contradiction(2, {0, 1, 1}, MAXIMIZE);
    contradiction.add_constraint({0, 1, 1}, 10, LESSEQ);
    contradiction.add_constraint({0, 1, 1}, 20, GREATEQ);

    assert(contradiction.solve(result) == INFEASIBLE);

    // Header example, both LPs built before either is solved
    Simplex header(2, {0, 3, 2}, MAXIMIZE);
    Simplex diet(2, {0, 1, 2}, MINIMIZE);
    header.add_constraint({0, 1, 1}, 4, LESSEQ);
    diet.add_constraint({0, 1, 1}, 3, GREATEQ);
    header.add_constraint({0, 2, 1}, 5, LESSEQ);
    diet.add_constraint({0, 1, -1}, 1, EQUAL);

    assert(header.solve(result) == FEASIBLE);
    assert(abs(result - 9) < 1e-9);
    assert(abs(header.val[1] - 1) < 1e-9 && abs(header.val[2] - 3) < 1e-9);

    assert(diet.solve(result) == FEASIBLE);
    assert(abs(result - 4) < 1e-9);
    assert(abs(diet.val[1] - 2) < 1e-9 && abs(diet.val[2] - 1) < 1e-9);

    // Unbounded, and no constraints at all
    Simplex runaway(2, {0, 1, 0}, MAXIMIZE);
    runaway.add_constraint({0, 1, -1}, 1, LESSEQ);
    assert(runaway.solve(result) == UNBOUNDED);

    Simplex free_lunch(2, {0, 1, 0}, MAXIMIZE);
    assert(free_lunch.solve(result) == UNBOUNDED);

    Simplex idle(2, {0, -1, -3}, MAXIMIZE);
    assert(idle.solve(result) == FEASIBLE);
    assert(abs(result) < 1e-9 && abs(idle.val[1]) < 1e-9 && abs(idle.val[2]) < 1e-9);

    return 0;
}
