/***
 *
 * Simplex
 * Maximizes or minimizes a linear objective subject to linear constraints on variables x_j >= 0
 *
 * Complexity: O(n m) time per pivot and O(n m) memory, the pivot count is exponential in the worst case
 *
 * Simplex lp(objective, sense): an LP in n = objective.size() variables, sense is Simplex::MAXIMIZE or Simplex::MINIMIZE
 * BasicSimplex<long double> lp(objective, sense): the same in long double, 2-3x slower, resolving 1e-13 instead of 1e-9
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
 * Each pivot picks steepest edge while pricing it costs less than the pivot itself, otherwise the largest reduced cost
 * Rows are scaled internally, so mixing very large and very small coefficients is fine
 * The answer is unreliable only when moving a limit by about 1e-7 of its row changes it, as with nearly parallel rows
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Outside the template so Simplex::MAXIMIZE also serves BasicSimplex<long double>
struct SimplexOptions{
    enum Cmp{LESSEQ, EQUAL, GREATEQ};
    enum Sense{MINIMIZE = -1, MAXIMIZE = 1};
    enum Status{INFEASIBLE, FEASIBLE, UNBOUNDED};
};

template <typename Float>
struct BasicSimplex : SimplexOptions{
    static constexpr bool WIDE = is_same<Float, long double>::value;
    /// EPS decides pivots and feasibility, DROP only prunes elimination: pruning at EPS broke rows once variables grew large
    static constexpr Float EPS = WIDE ? 1e-13L : 1e-9, DROP = WIDE ? 1e-16L : 1e-12;

    Float value = 0;
    vector<Float> val;
    int n, m, pivot_width;
    Sense sense;
    vector<vector<Float>> ar;
    vector<Float> rhs;
    vector<int> nonzero, nonbasic, basic;

    BasicSimplex(const vector<Float>& objective, Sense sense)
        : val(objective.size(), 0), n(objective.size()), m(0), pivot_width(n + 1), sense(sense), ar(1, vector<Float>(n + 1, 0)), rhs(1, 0), nonzero(n + 1), nonbasic(n + 1){
        for (int j = 0; j < n; j++) ar[0][j + 1] = objective[j] * sense;
    }

    void add_constraint(const vector<Float>& coef, Float lim, Cmp cmp){
        assert((int)coef.size() == n);
        if (cmp != GREATEQ) add_row(coef, lim, 1);
        if (cmp != LESSEQ) add_row(coef, lim, -1);
    }

    Status solve(){
        assert(basic.empty());  /// the tableau is pivoted in place, so a second call would read garbage

        /// Variable ids: 0 is the phase 1 artificial, 1..n the inputs, n + i the slack of row i
        basic.assign(m + 1, 0);
        for (int j = 0; j <= n; j++) nonbasic[j] = j;
        for (int i = 1; i <= m; i++) basic[i] = i + n;

        /// Phase 1: an artificial x0 in column 0 relaxes every row, maximizing -x0 reaches a feasible basis
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

            /// x0 can stay basic within EPS and grow in phase 2, so pivot it out; a row with no pivot is redundant
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
        Float scale = 0;
        for (Float c : coef) scale = max(scale, abs(c));
        if (scale == 0) scale = 1;
        m++;
        rhs.push_back(lim * sign / scale);
        ar.emplace_back(n + 1, 0);
        for (int j = 0; j < n; j++) ar[m][j + 1] = coef[j] * sign / scale;
    }

    bool optimize(int z){
        /// Bland's rule after a degenerate pivot, since only degenerate pivots can cycle
        for (bool bland = false; ; ){
            int x = 0, y = -1;
            Float v = 0;
            if (bland){
                for (int j = 0; j <= n; j++){
                    if (ar[z][j] > EPS && (y < 0 || nonbasic[j] < nonbasic[y])) y = j;
                }
            }
            else{
                /// Steepest edge costs O(m) per candidate, so it runs only while that is cheaper than the last pivot, else Dantzig
                int candidates = 0;
                for (int j = 0; j <= n; j++) candidates += ar[z][j] > EPS;
                bool steepest = candidates <= pivot_width;
                Float best = 0;
                for (int j = 0; j <= n; j++){
                    if (ar[z][j] <= EPS) continue;
                    Float len = 1;
                    for (int i = 1; i <= m && steepest; i++) len += ar[i][j] * ar[i][j];
                    if (ar[z][j] * ar[z][j] / len > best) best = ar[z][j] * ar[z][j] / len, y = j;
                }
            }
            if (y < 0) return true;

            for (int i = 1; i <= m; i++){
                if (ar[i][y] <= EPS) continue;
                Float u = rhs[i] / ar[i][y];
                /// Near ties take the largest pivot to limit rounding, or the lowest id under Bland's rule
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
        pivot_width = len;

        /// Row 0 is the objective and row m + 1 the phase 1 objective while it exists, each rhs holds minus its value
        for (int i = 0; i < (int)ar.size(); i++){
            vector<Float>& cur = ar[i];
            if (i == x || abs(cur[y]) <= DROP) continue;
            rhs[i] -= cur[y] * rhs[x], v = cur[y], cur[y] = 0;
            for (int k = 0; k < len; k++) cur[nonzero[k]] -= v * row[nonzero[k]];
        }
    }
};

using Simplex = BasicSimplex<double>;

/// Groups 0..6 are a, b, c, ab, ac, bc, abc and variable 7 is the hour count T
const int GROUPS = 7, VARS = 8;
const int USES[4][GROUPS] = {{1, 0, 0, 1, 1, 0, 1}, {0, 1, 0, 1, 0, 1, 1}, {0, 0, 1, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1, 1}};

int best;
long long need[GROUPS];
vector<int> order;
vector<double> lo(VARS, 0), hi(VARS, -1);

void branch_and_bound(){
    vector<double> objective(VARS, 0);
    objective[GROUPS] = 1;
    Simplex lp(objective, Simplex::MINIMIZE);

    /// Musketeer k fights at most once per hour, and pairs or the trio exclude each other within an hour
    for (int k = 0; k < 4; k++){
        vector<double> coef(USES[k], USES[k] + GROUPS);
        coef.push_back(-1);
        lp.add_constraint(coef, 0, Simplex::LESSEQ);
    }

    /// Hall's condition on strengths: the criminals needing at least order[i] get enough fights from that group or stronger
    for (int i = 0; i < GROUPS; i++){
        vector<double> coef(VARS, 0);
        long long criminals = 0;
        for (int j = i; j < GROUPS; j++) coef[order[j]] = 1, criminals += need[j];
        lp.add_constraint(coef, criminals, Simplex::GREATEQ);
    }

    for (int j = 0; j < VARS; j++){
        vector<double> coef(VARS, 0);
        coef[j] = 1;
        if (lo[j] > 0) lp.add_constraint(coef, lo[j], Simplex::GREATEQ);
        if (hi[j] >= 0) lp.add_constraint(coef, hi[j], Simplex::LESSEQ);
    }

    if (lp.solve() != Simplex::FEASIBLE || ceil(lp.value - 1e-6) >= best) return;
    for (int j = 0; j < VARS; j++){
        double down = floor(lp.val[j] + 1e-6);
        if (lp.val[j] - down < 1e-6) continue;

        double old_lo = lo[j], old_hi = hi[j];
        lo[j] = down + 1;
        branch_and_bound();
        lo[j] = old_lo, hi[j] = down;
        branch_and_bound();
        hi[j] = old_hi;
        return;
    }
    best = llround(lp.value);
}

int main(){
    int n;
    long long a, b, c;
    scanf("%d %lld %lld %lld", &n, &a, &b, &c);
    long long strength[GROUPS], values[GROUPS] = {a, b, c, a + b, a + c, b + c, a + b + c};
    for (int i = 0; i < GROUPS; i++) order.push_back(i);
    sort(order.begin(), order.end(), [&](int x, int y){
        return make_pair(values[x], x) < make_pair(values[y], y);
    });
    for (int i = 0; i < GROUPS; i++) strength[i] = values[order[i]];

    for (int i = 0; i < n; i++){
        long long t;
        scanf("%lld", &t);
        if (t > a + b + c){
            puts("-1");
            return 0;
        }
        need[lower_bound(strength, strength + GROUPS, t) - strength]++;
    }

    best = n;
    branch_and_bound();
    printf("%d\n", best);
    return 0;
}
