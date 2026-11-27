#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/simplex.cpp"
#undef main

typedef __int128 i128;
using Float = double;

struct Constraint{
    vector<long long> a;
    long long b;
    Simplex::Cmp cmp;
};

struct Fraction{
    i128 num, den;  /// den > 0
    bool operator<(const Fraction& o) const { return num * o.den < o.num * den; }
    long double value() const { return (long double)num / (long double)den; }
};

i128 det(vector<vector<i128>> m){
    int n = m.size();
    i128 res = 1;
    /// Bareiss keeps every intermediate an exact integer, entries here stay far below the 128-bit range
    i128 prev = 1;

    for (int k = 0; k < n; k++){
        int p = k;
        while (p < n && m[p][k] == 0) p++;
        if (p == n) return 0;
        if (p != k) swap(m[p], m[k]), res = -res;
        for (int i = k + 1; i < n; i++){
            for (int j = k + 1; j < n; j++) m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) / prev;
        }
        prev = m[k][k];
    }

    return res * m[n - 1][n - 1];
}

/// Best objective over vertices of {constraints, 0 <= x <= box}, or nothing when that region is empty
bool brute(int n, const vector<long long>& c, const vector<Constraint>& cons, Simplex::Sense sense, long long box, Fraction& best){
    vector<pair<vector<long long>, long long>> planes;
    for (auto& k : cons) planes.push_back({k.a, k.b});
    for (int i = 0; i < n; i++){
        vector<long long> e(n, 0);
        e[i] = 1;
        planes.push_back({e, 0}), planes.push_back({e, box});
    }

    bool found = false;
    int total = planes.size();
    vector<int> pick(n);
    function<void(int, int)> choose = [&](int start, int depth){
        if (depth == n){
            vector<vector<i128>> m(n, vector<i128>(n));
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) m[i][j] = planes[pick[i]].first[j];
            i128 d = det(m);
            if (d == 0) return;
            vector<i128> num(n);
            for (int j = 0; j < n; j++){
                auto mj = m;
                for (int i = 0; i < n; i++) mj[i][j] = planes[pick[i]].second;
                num[j] = det(mj);
            }
            if (d < 0){
                d = -d;
                for (auto& v : num) v = -v;
            }

            for (int j = 0; j < n; j++) if (num[j] < 0 || num[j] > box * d) return;
            for (auto& k : cons){
                i128 s = 0;
                for (int j = 0; j < n; j++) s += k.a[j] * num[j];
                i128 rhs = (i128)k.b * d;
                if ((k.cmp == Simplex::LESSEQ && s > rhs) || (k.cmp == Simplex::GREATEQ && s < rhs) || (k.cmp == Simplex::EQUAL && s != rhs)) return;
            }

            Fraction value{0, d};
            for (int j = 0; j < n; j++) value.num += c[j] * num[j];
            if (sense == Simplex::MINIMIZE) value.num = -value.num;
            if (!found || best < value) best = value;
            found = true;
            return;
        }
        for (int i = start; i < total; i++){
            pick[depth] = i;
            choose(i + 1, depth + 1);
        }
    };

    choose(0, 0);
    if (found && sense == Simplex::MINIMIZE) best.num = -best.num;
    return found;
}

/// Phase 1 leaves an artificial basic at 2e-10, inside EPS, which grows in phase 2 and relaxes x0 >= 1000 away
void artificial_left_basic(){
    Simplex lp({1}, Simplex::MINIMIZE);
    lp.add_constraint({0.001}, 1, Simplex::GREATEQ);
    lp.add_constraint({1}, 1000 - 2e-7, Simplex::LESSEQ);

    int status = lp.solve();
    assert(status == Simplex::INFEASIBLE || (status == Simplex::FEASIBLE && abs(lp.value - 1000) < 1e-6));
}

/// Degenerate LPs that cycle forever without anti-cycling, the alarm turns a hang into a failure
void degenerate_cycling(){
    alarm(2);

    /// Beale (1955): Dantzig's rule with first-row ratio ties cycles here, optimum 1/20 at x0 = 1/25, x2 = 1
    Simplex beale({0.75, -150, 0.02, -6}, Simplex::MAXIMIZE);
    beale.add_constraint({0.25, -60, -0.04, 9}, 0, Simplex::LESSEQ);
    beale.add_constraint({0.5, -90, -0.02, 3}, 0, Simplex::LESSEQ);
    beale.add_constraint({0, 0, 1, 0}, 1, Simplex::LESSEQ);
    assert(beale.solve() == Simplex::FEASIBLE);
    assert(abs(beale.value - 0.05) < 1e-9);
    assert(abs(beale.val[0] - 0.04) < 1e-9 && abs(beale.val[2] - 1) < 1e-9);

    /// Each cycles under a different broken tie rule, rows are {coefficients..., rhs} plus sum(x) <= 10
    vector<tuple<vector<Float>, vector<vector<Float>>, Float>> cases = {
        {{3, -2, 2, 6, 6}, {{-2, 0, 0, -1, -3, 0}, {0, 3, 0, 0, 0, 0}, {0, 3, 4, 0, -3, 1}, {0, -4, 1, 4, 2, 0}, {0, -2, -1, -2, 4, 1}, {1, 0, 0, 1, 3, 0}, {4, 2, 2, 0, 2, 0}}, 0},
        {{3, -1, -1, -1, -4, 0, -2}, {{0, -4, 0, 0, 2, 1, 2, 0}, {4, 3, 0, -1, 0, -1, 1, 0}, {2, 0, -2, -4, 2, 0, 2, 0}, {-2, 0, 0, 4, 3, -2, 0, 0}, {0, 4, 1, 3, 3, -2, 2, 0}, {1, 0, 2, 0, 0, -4, 3, 0}, {-3, -4, -3, -3, 0, 0, 0, 0}}, 0},
        {{-1, 5, -3, -1, -1, 3, 2}, {{-2, 0, -3, 1, 1, 2, 1, 3}, {0, 2, 2, -2, 0, -2, 0, 0}, {3, 4, -2, 4, 0, 4, -3, 0}, {2, 1, -3, 0, 3, 0, -4, 0}, {0, -2, 4, 2, 0, 0, 0, 0}, {-3, 0, 1, 2, 1, 0, 0, 0}, {0, -1, -3, 0, 0, -4, -3, 0}}, 123.0 / 7},
        {{-3, 2, -1, -2, -3, 3, 4}, {{-4, -4, 0, -1, -2, -1, 4, 0}, {-1, 2, 1, 0, -2, 2, -2, 0}, {0, 2, 3, -1, 0, 0, -1, 0}, {4, 0, -2, -3, -3, 0, 3, 0}, {2, -3, 0, 1, -2, 1, 4, 0}}, 65.0 / 8},
        /// Cycles when degenerate ratio ties take the largest pivot instead of Bland's lowest id
        {{-3, 5, 1, 3, 3, 4, 3}, {{1, -2, -3, 2, 2, 0, 0, 0}, {0, 0, -1, -2, 3, 0, 3, 0}, {-4, 1, -2, 3, -3, -1, -4, 0}, {-2, 0, -2, 0, -2, -2, 3, 0}, {0, 4, 2, 0, 0, 0, -1, 0}}, 40},
    };
    for (auto& [obj, rows, opt] : cases){
        int n = obj.size();
        Simplex lp(obj, Simplex::MAXIMIZE);
        for (auto& r : rows) lp.add_constraint(vector<Float>(r.begin(), r.begin() + n), r[n], Simplex::LESSEQ);
        lp.add_constraint(vector<Float>(n, 1), 10, Simplex::LESSEQ);
        assert(lp.solve() == Simplex::FEASIBLE);
        assert(abs(lp.value - opt) < 1e-9);
    }

    alarm(0);
}

/// Primal max c.x, Ax <= b and dual min b.y, A^T y >= c as two live instances, strong duality pins statuses and optima
void duality(long long iterations, int lo, int hi, bool planted){
    for (long long it = 0; it < iterations; it++){
        int n = stress::rand_int(lo, hi), m = stress::rand_int(lo, hi);
        vector<vector<long long>> a(m, vector<long long>(n));
        vector<long long> b(m), c(n), x0(n, 0), y0(m, 0);
        for (auto& r : a) for (auto& x : r) x = stress::rand_int(0, 3) ? stress::rand_int(-4, 6) : 0;
        /// b, c from x0, y0 >= 0 keep both sides feasible and bounded; arbitrary LPs add infeasible and unbounded pairs
        if (planted){
            for (auto& x : x0) x = stress::rand_int(0, 1) ? stress::rand_int(0, 5) : 0;
            for (auto& y : y0) y = stress::rand_int(0, 1) ? stress::rand_int(0, 5) : 0;
        }
        for (int i = 0; i < m; i++){
            b[i] = planted ? (stress::rand_int(0, 1) ? stress::rand_int(0, 10) : 0) : stress::rand_int(-10, n + 10);
            for (int j = 0; j < n; j++) b[i] += a[i][j] * x0[j];
        }
        for (int j = 0; j < n; j++){
            c[j] = planted ? (stress::rand_int(0, 1) ? -stress::rand_int(0, 10) : 0) : stress::rand_int(-m - 10, 10);
            for (int i = 0; i < m; i++) c[j] += a[i][j] * y0[i];
        }

        vector<Float> row(n), col(m);
        Simplex primal(vector<Float>(c.begin(), c.end()), Simplex::MAXIMIZE), dual(vector<Float>(b.begin(), b.end()), Simplex::MINIMIZE);
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++) row[j] = a[i][j];
            primal.add_constraint(row, b[i], Simplex::LESSEQ);
        }
        for (int j = 0; j < n; j++){
            for (int i = 0; i < m; i++) col[i] = a[i][j];
            dual.add_constraint(col, c[j], Simplex::GREATEQ);
        }

        alarm(10);  /// past n, m = 50 a phase 1 without a ratio test hangs, the alarm turns that into a failure
        int primal_status = primal.solve(), dual_status = dual.solve();
        alarm(0);
        Float primal_value = primal.value, dual_value = dual.value;
        if (planted) assert(primal_status == Simplex::FEASIBLE);
        assert((primal_status == Simplex::FEASIBLE) == (dual_status == Simplex::FEASIBLE));
        assert(primal_status != Simplex::UNBOUNDED || dual_status == Simplex::INFEASIBLE);
        assert(dual_status != Simplex::UNBOUNDED || primal_status == Simplex::INFEASIBLE);
        if (primal_status != Simplex::FEASIBLE) continue;

        long double tol = 1e-6L * max(1.0L, fabsl(primal_value));
        assert(fabsl(primal_value - dual_value) <= tol);

        long double objective = 0;
        for (int j = 0; j < n; j++){
            assert(primal.val[j] >= -1e-7);
            objective += c[j] * primal.val[j];
        }
        assert(fabsl(objective - primal_value) <= tol);
        for (int i = 0; i < m; i++){
            long double s = 0;
            for (int j = 0; j < n; j++) s += a[i][j] * primal.val[j];
            assert(s <= b[i] + 1e-6);
        }

        objective = 0;
        for (int i = 0; i < m; i++){
            assert(dual.val[i] >= -1e-7);
            objective += b[i] * dual.val[i];
        }
        assert(fabsl(objective - dual_value) <= tol);
        for (int j = 0; j < n; j++){
            long double s = 0;
            for (int i = 0; i < m; i++) s += a[i][j] * dual.val[i];
            assert(s >= c[j] - 1e-6);
        }
    }
}

/// UVa 1298: nearly parallel rows of 1e5 coefficients, infeasible, but FEASIBLE to a lowest-id rule among near ties
void uva_1298_rows(){
    int s[7][3] = {{2470, 6466, 793}, {2472, 6468, 794}, {2470, 6472, 791}, {2472, 6469, 789}, {9991, 9998, 9999}, {10000, 9999, 10000},
                   {9994, 9996, 10000}};
    Simplex lp({0, 0, 0}, Simplex::MAXIMIZE);
    for (int j = 0; j < 6; j++){
        vector<Float> c(3);
        for (int k = 0; k < 3; k++) c[k] = 1e8 / s[j][k] - 1e8 / s[6][k];
        lp.add_constraint(c, 1, Simplex::GREATEQ);
    }
    assert(lp.solve() == Simplex::INFEASIBLE);
}

/// Rows of very different size: unscaled, the optimum value stays right while y breaks 600000y >= 1100000
void mixed_magnitude_rows(){
    Simplex lp({-4, -3}, Simplex::MAXIMIZE);
    lp.add_constraint({0, 0}, 800, Simplex::LESSEQ);
    lp.add_constraint({0, 600000}, 1100000, Simplex::GREATEQ);
    lp.add_constraint({-40, 0}, -120, Simplex::LESSEQ);
    lp.add_constraint({500000, 0}, 300000, Simplex::GREATEQ);
    assert(lp.solve() == Simplex::FEASIBLE);
    assert(abs(lp.value + 17.5) < 1e-9);
    assert(abs(lp.val[0] - 3) < 1e-9);
    assert(abs(lp.val[1] - 11.0 / 6) < 1e-9);
}

/// Nearly parallel rows, infeasible even with limits moved by 1e-7 of their row, FEASIBLE to a lowest-id near-tie rule
void nearly_parallel_rows(){
    Simplex lp({-4, 4, -2}, Simplex::MINIMIZE);
    lp.add_constraint({0, -1, -2}, 3, Simplex::LESSEQ);
    lp.add_constraint({89997, -30003, -90002}, 30000, Simplex::GREATEQ);
    lp.add_constraint({89999, -29997, -90001}, 30000, Simplex::EQUAL);
    lp.add_constraint({90002, -30001, -90002}, 4, Simplex::GREATEQ);
    assert(lp.solve() == Simplex::INFEASIBLE);
}

/// UVa 12512: grid cross sums as opposite row pairs, where a 2e-13 pivot on a Bland tie broke rows or never ended
template <typename F>
void uva_12512_equalities(){
    struct Case{ int c, d, f, target; vector<int> g; };
    vector<Case> cases = {
        {51, 63, 9, 66972, {5, 1, 3, 9, 8, 3, -1, 5, 4, 3, 3, -7, 8, -7, 6, 0, 6, 6, 3, 3, -10, 0, -9, -1, 9, 3, 2, 9, 0, -1, 8, -8, -3, -7, 4, -3,
                            6, 2, -9, 5, -7, 6, -3, 7, -1, 2, -3, -8, 1, 7, 3, 8, -5, -3, 10, 3, 8, -7, -8, -7, 6, -6, -7, 7, 8, 4, -10, -1, -3, -1,
                            -10, 0, -5, -8, 4, 4, 6, 2, -6, 1, 8, -3, -9, -1, -10, 4, 0, -1, 0, 5, -6, 5, -4, -4, -1, -9, -4, 0, 1, 5, -10, -5, 10, 7,
                            -4, 7, 10, -3, -4, -9, 10, -7, -1, 0, -2, 3, -6, 0, 4, 3, -8, -5, -3, -1, -7, -4, -2, -5, -3, 6, -1, -2, 1, 6, -7, 8, 6, -3,
                            -1, -5, -3, -2, -6, 10}},
        {85, 11, 4, 60904, {4, 6, -8, 6, -9, -10, 9, 8, -1, 7, -2, 3, 8, -8, -5, -1, -4, -3, -5, -2, 7, -10, -5, -6, -9, 9, 4, 8, 2, -10, 5, -9, -5, -8,
                            7, 10, 1, 0, 5, -9, 6, -1, 6, 10, 6, 7, -5, 9, -5, -7, 9, 3, -3, 4, 7, -2, -3, -2, 5, 7, -3, 8, 9, 5, -2, -9, -4, 9, -9, 8,
                            4, -5, 0, -6, 9, -10, 9, 5, 3, 6, 4, -2, -10, -8, -7, -5, 8, 10, 7, 8, -3, 9, -6, -9, -1, 3, -9, 0, 9, 6, -2, 8, 5, -8, 10,
                            6, -5, 2, -2, -3, -5, 5, 9, 4, -9, 0, 3, 4, 1, -4, -9, -6, 5, 9, -6, -6, 4, 6, -9, -10, -3, 3, 6, -6, -8, 2, 7, -1, 5, 3, -8,
                            -9, 9, -3}},
    };

    alarm(10);
    const int n = 12;
    for (auto& t : cases){
        vector<F> cost(2 * n * n);
        for (int i = 0; i < n * n; i++) cost[2 * i] = t.c, cost[2 * i + 1] = t.d;

        BasicSimplex<F> lp(cost, Simplex::MINIMIZE);
        vector<vector<F>> rows;
        vector<F> lim;
        for (int x = 0; x < n; x++){
            for (int y = 0; y < n; y++){
                vector<F> coef(2 * n * n, 0);
                long long sum = 0;
                for (int i = 0; i < n; i++){
                    for (int j = 0; j < n; j++){
                        if (i != x && j != y) continue;
                        sum += t.g[i * n + j], coef[2 * (i * n + j)] = 1, coef[2 * (i * n + j) + 1] = -1;
                    }
                }
                lp.add_constraint(coef, (2 * n - 1) * t.f - sum, Simplex::LESSEQ);
                lp.add_constraint(coef, (2 * n - 1) * t.f - sum, Simplex::GREATEQ);
                rows.push_back(coef), lim.push_back((2 * n - 1) * t.f - sum);
            }
        }

        assert(lp.solve() == Simplex::FEASIBLE);
        assert(abs(lp.value - t.target) < 1e-6 * t.target);
        for (int r = 0; r < (int)rows.size(); r++){
            long double s = 0;
            for (int j = 0; j < 2 * n * n; j++) s += rows[r][j] * lp.val[j];
            assert(abs(s - lim[r]) < 1e-6);
        }
    }
    alarm(0);
}

/// Roundoff degenerate limits must count as exactly 0, or ties break by noise instead of Bland's lowest id and cycle
void roundoff_degenerate_rows(){
    alarm(2);
    vector<Float> objective = {57, 36, 4, -18, -16, 45, 96, -1, 35, 59, 63, 100, 57};
    /// {coefficients..., limit}, optimum 562 by HiGHS
    vector<vector<Float>> rows = {
        {-2, 0, 0, 0, 5, 0, 5, 0, 1, -3, 6, 6, 4, 22},
        {5, 0, -4, 6, 1, 1, 0, 3, -2, 4, 0, 5, -4, -4},
        {-4, -2, 0, 5, 0, -2, -4, 3, 0, 3, 3, 0, 6, -12},
        {6, -4, 3, -3, 2, 0, 0, 2, 4, 4, 2, -2, -4, -4},
        {3, 4, -3, 2, 2, 0, 6, -2, 1, 3, 0, -3, 3, 38},
        {-3, -4, -2, 3, 0, 0, 0, -1, 3, 6, 5, -1, 0, -6},
        {6, 6, -4, -4, -1, -1, -1, -2, 4, 0, 0, 4, 0, 22},
        {0, 4, 0, 0, -2, -1, -4, 0, -3, 5, 4, 5, -3, -10},
        {4, 0, 0, -2, 1, 6, 0, 3, -1, -4, -2, 5, 4, -2},
        {0, 3, 0, 0, 0, 5, 5, 4, 0, 0, 2, 5, 2, 29},
        {-1, -2, 2, 1, -3, -1, 0, 0, -1, -3, 0, -2, 5, -8}
    };
    int n = objective.size();
    Simplex lp(objective, Simplex::MAXIMIZE);
    for (auto& r : rows) lp.add_constraint(vector<Float>(r.begin(), r.begin() + n), r[n], Simplex::LESSEQ);
    assert(lp.solve() == Simplex::FEASIBLE);
    assert(abs(lp.value - 562) < 1e-6);
    alarm(0);
}

/// Codeforces 605C in dual form has corners near 1e-6, where an absolute 1e-9 ratio tolerance ties ratios 1e-12 apart
void tiny_ratio_rows(long long iterations){
    for (long long it = 0; it < iterations; it++){
        int n = stress::rand_int(2, 60);
        long long lo = it % 2 ? 900000 : 1, p = stress::rand_int(1, 1000000), q = stress::rand_int(1, 1000000);
        vector<long long> a(n), b(n);
        for (int i = 0; i < n; i++) a[i] = stress::rand_int(lo, 1000000), b[i] = stress::rand_int(lo, 1000000);

        /// The primal, min sum t_i with sum a_i t_i >= p and sum b_i t_i >= q, is optimal on one project or a pair
        long double best = 1e30L;
        for (int i = 0; i < n; i++){
            best = min(best, max((long double)p / a[i], (long double)q / b[i]));
            for (int j = 0; j < i; j++){
                long double d = (long double)a[i] * b[j] - (long double)a[j] * b[i];
                if (d == 0) continue;
                long double ti = ((long double)p * b[j] - (long double)q * a[j]) / d, tj = ((long double)q * a[i] - (long double)p * b[i]) / d;
                if (ti >= 0 && tj >= 0) best = min(best, ti + tj);
            }
        }

        Simplex lp({(Float)p, (Float)q}, Simplex::MAXIMIZE);
        for (int i = 0; i < n; i++) lp.add_constraint({(Float)a[i], (Float)b[i]}, 1, Simplex::LESSEQ);
        assert(lp.solve() == Simplex::FEASIBLE);
        assert(fabsl(lp.value - best) <= 1e-9L * best);
    }
}

/// The tableau is pivoted in place, so a second solve() on the same instance must abort instead of answering
void solve_twice_aborts(){
    stress::expect_assert("linear_algebra/simplex.cpp", []{
        Simplex encore({3, 2}, Simplex::MAXIMIZE);
        encore.add_constraint({1, 1}, 4, Simplex::LESSEQ);
        encore.solve();
        encore.solve();
    });
}

/// After INFEASIBLE there is no feasible corner to start from, so resolve must abort, as must resolve before solve
void resolve_without_basis_aborts(){
    stress::expect_assert("linear_algebra/simplex.cpp", []{
        Simplex lp({1, 1}, Simplex::MAXIMIZE);
        lp.add_constraint({1, 1}, 1, Simplex::LESSEQ);
        lp.add_constraint({1, 1}, 2, Simplex::GREATEQ);
        lp.solve();
        lp.resolve({1, 0});
    });
    stress::expect_assert("linear_algebra/simplex.cpp", []{
        Simplex lp({1, 1}, Simplex::MAXIMIZE);
        lp.add_constraint({1, 1}, 1, Simplex::LESSEQ);
        lp.resolve({1, 0});
    });
}

/// Larger LPs: resolve over one tableau must match a fresh solve of the same constraints with each new objective
template <typename Float>
void resolve_matches_fresh(long long iterations){
    for (long long it = 0; it < iterations; it++){
        int n = stress::rand_int(5, 40), m = stress::rand_int(1, 40);
        auto sense = stress::rand_int(0, 1) ? Simplex::MAXIMIZE : Simplex::MINIMIZE;
        vector<vector<Float>> rows(m, vector<Float>(n));
        vector<Float> lims(m);
        vector<Simplex::Cmp> cmps(m);
        for (int i = 0; i < m; i++){
            for (auto& x : rows[i]) x = stress::rand_int(0, 3) ? stress::rand_int(-9, 9) : 0;
            lims[i] = stress::rand_int(-30, 100);
            int t = stress::rand_int(0, 7);
            cmps[i] = t < 5 ? Simplex::LESSEQ : t < 7 ? Simplex::GREATEQ : Simplex::EQUAL;
        }
        auto build = [&](const vector<Float>& c){
            BasicSimplex<Float> lp(c, sense);
            for (int i = 0; i < m; i++) lp.add_constraint(rows[i], lims[i], cmps[i]);
            return lp;
        };
        auto objective = [&](){
            vector<Float> c(n);
            for (auto& x : c) x = stress::rand_int(-9, 9);
            return c;
        };

        BasicSimplex<Float> lp = build(objective());
        if (lp.solve() == Simplex::INFEASIBLE) continue;
        for (int q = 0; q < 5; q++){
            vector<Float> c = objective();
            BasicSimplex<Float> fresh = build(c);
            int want = fresh.solve();
            assert(lp.resolve(c) == want);
            if (want == Simplex::FEASIBLE) assert(fabsl(lp.value - fresh.value) <= 1e-6L * max((long double)1, fabsl(fresh.value)));
        }
    }
}

/// Nightly seed 45993096780: a skipped 1.6e-10 binding entry left the dual infeasible, 1.3e-5 below the optimum
void nightly_skipped_pivot(){
    int n = 32, m = 2 * n;
    const int cost[] = {36, 14, 40, 75, 46, 20, 18, 23, 78, 46, 51, 82, 5, 56, 28, 96, 71, 66, 85, 78, 34, 83, 67, 38, 51, 91, 97, 15, 58, 58, 54, 32};
    const int col_exp[] = {-1, 3, -3, -3, -4, -6, 1, -7, -1, 6, -3, 9, 4, -2, 5, -3, -6, 10, 10, -4, 0, 9, -2, 2, -3, 6, -1, 5, -7, 8, 3, -6};
    const int row_exp[] = {-2, 10, -6, -2, -12, 4, 8, 3, 10, 3, -6, 7, -7, -3, -11, 6, -13, 13, 7, 0, 2, -13, 10, -13, -13, -4, 5, 1, 2, 10, 3, -13};
    const int lim[] = {1572, 133, 6246, 16, 8518, 44, 5126, 442, 4465, 355, 2811, 59, 9414, 185, 5102, 447, 5271, 194, 7145, 97, 5177, 246, 3133, 109, 1376, 82, 5893, 307, 7090, 486, 2178, 10};
    const int coef[32][32] = {
        {13, 70, 67, 28, -40, -27, -35, 64, 57, 2, 31, 51, 88, -30, -46, -34, 10, -9, 82, 88, 65, 31, 98, -16, 35, -31, 4, 85, 67, -3, 72, 51},
        {51, 40, 69, -6, 38, 98, -31, 17, 86, -37, -31, 68, -12, 42, -8, 82, 52, 8, 93, -21, -17, 57, 75, -43, -26, -7, 5, 28, 2, 13, 52, 69},
        {-8, 55, -13, -50, 48, 97, -27, 0, 76, 28, -49, 3, -38, -11, -25, 96, 36, -39, 78, 11, -25, 55, 18, 35, 30, 43, 61, 6, -22, 98, 46, -6},
        {-10, -50, 58, -19, 34, 36, 90, 72, -38, 37, 63, -8, -13, 1, -35, 88, 24, 55, 86, -44, 98, 35, 98, 66, -31, -5, 84, -25, 29, 12, -2, 62},
        {2, 75, 86, 47, 45, -47, 29, 55, 93, 9, 6, 33, -1, 81, 61, 21, -41, 52, -39, 65, 96, 16, -4, 25, 63, -3, 23, -36, 61, -26, 100, -37},
        {15, -14, -22, 44, 36, 91, 7, -29, -44, 63, 54, -12, 48, -50, 26, -40, 99, 58, 53, 48, -45, 27, -48, -34, 72, 82, -19, 50, 97, 62, 13, 17},
        {39, 80, 21, 60, 37, 33, 61, -40, 38, 72, 4, 50, 74, -12, -11, -32, 60, 44, 44, 66, -7, -46, 21, 57, -1, -32, 80, 76, 85, -40, -15, -42},
        {51, 1, -17, 50, -44, 95, 61, 49, 17, 26, 63, 19, 95, 19, -14, 6, -30, 24, 36, 61, 34, -22, -4, 21, 70, 16, -47, 95, -4, 28, 9, 12},
        {79, -40, -44, 13, 54, -34, 17, 43, 37, 43, -6, -24, -13, 72, 77, -48, 45, -20, 33, 78, 38, 73, 27, 85, 49, 61, -31, 71, 40, 53, 65, 50},
        {68, 88, -15, -1, -43, 20, 89, -24, 9, 97, 33, -2, 22, -24, 43, 5, 15, -40, 92, 58, 71, 85, 90, -2, 58, 1, 8, 64, 60, -45, -28, -25},
        {-43, 95, 9, 62, 64, 32, -38, -23, 86, 90, -30, -28, -22, 14, 51, 0, 22, 6, 93, 25, 88, -40, 69, -21, 55, 5, 4, -20, 15, 91, -21, 1},
        {49, 20, 27, -6, 40, 97, 37, 81, 56, 61, 21, 87, 23, -34, 40, 95, -8, 42, -38, 63, 46, -45, 23, 58, 32, -31, 44, -50, 95, 49, -45, 95},
        {84, 22, 73, 78, -22, -33, -26, -9, -47, -45, 6, 88, -45, -41, 9, 52, 50, 7, -29, -30, -17, -27, 44, -49, 68, 86, -21, 31, 38, 44, 73, 43},
        {16, 14, 6, 26, 33, 41, -10, 45, 35, -17, 80, -23, 88, 57, 60, -9, 7, -23, 88, 22, 58, 89, 82, 28, 46, 97, 2, 16, -18, 19, -20, 77},
        {-47, 18, -43, 51, 17, 34, 72, 95, 44, 88, -28, -20, 33, -40, 20, -30, 27, -1, 15, -42, -36, -19, 44, 89, 65, -4, 66, 7, -1, 84, 29, 13},
        {-12, 80, 51, 97, 93, 10, 62, -26, -21, 19, -24, -16, 6, 64, 97, 88, -33, 18, -14, 12, 32, 39, -35, -35, -41, 36, -6, -47, -45, 17, -24, -27},
        {4, 81, 29, 13, 49, 36, 10, 4, 37, 25, -27, -48, 53, 68, 54, 7, 92, -12, 95, 99, 14, -22, 83, 60, -3, 74, 6, 95, 42, -6, 23, 45},
        {-14, -26, 48, 24, 89, 87, -30, 85, 99, 34, 15, 41, 3, 68, 41, -44, 52, -20, -15, -14, -6, 99, -24, -21, -1, 3, -44, -11, 3, -11, 36, 49},
        {-21, -13, 44, -8, 8, -40, -19, -4, 37, 97, -35, -8, -48, -21, 20, 85, 73, 11, 78, 8, 70, 15, -33, 61, 14, 65, 12, -2, 79, 61, -43, 85},
        {-47, 49, -32, 67, 54, -30, 1, -4, -23, 41, 37, 30, -6, 70, 47, 27, -32, -48, 84, -21, 33, -18, -28, 8, -7, 4, -8, 42, 55, 77, 4, 85},
        {44, 84, 28, 50, 26, -16, 69, -20, 90, 90, 52, -37, -12, 42, -32, 53, 59, 14, -7, 52, -37, -12, 45, 2, -34, 73, 17, 77, 61, -13, -29, -37},
        {94, -8, -18, 12, 56, 26, 9, -33, 25, 37, 6, 0, -27, -12, -39, -12, 56, 31, 44, 23, -49, -33, 83, 97, 34, -31, 46, 90, -7, -35, 62, 41},
        {9, -35, -9, 17, 13, -8, -7, 98, 33, 2, -13, 48, 27, 88, -12, 75, 18, 80, -6, -22, -39, -22, 100, -33, -40, -49, 5, -27, 14, 12, 63, 31},
        {80, 75, -28, 28, 31, 63, 68, 93, 36, 88, -41, -2, 72, -29, 63, -37, -11, 71, 74, 83, 76, -21, 36, 43, -42, 10, -8, -35, 88, 96, 91, 31},
        {-45, -34, -8, -1, -27, -33, 1, 48, 16, 91, 96, 21, 50, -38, -44, -46, 41, -29, 15, 47, 52, 82, 35, -27, 14, 6, 44, 22, 60, 59, -7, 27},
        {-7, 96, 52, 97, 42, 58, -7, 16, 6, 19, 81, 56, -10, -4, -46, 11, -32, 79, 88, 46, 93, 24, 57, -40, -24, -46, 22, -43, 22, 70, -3, 79},
        {-47, 76, -6, 80, -48, 17, 2, 59, 71, -4, 46, 59, 1, -3, 91, 18, -36, -2, -16, 38, 70, 17, 94, -14, 49, 90, 96, 77, -22, 61, 58, 14},
        {46, -43, -23, 64, 11, -31, 4, 49, -40, 100, 14, 31, 80, 0, 58, 4, 84, 34, 72, -40, -6, -45, 2, 44, 15, 74, -43, -9, 19, -4, 0, -36},
        {-25, 3, 15, -7, 86, -39, 31, -7, 88, 59, 69, -30, 68, 26, 41, 84, 3, 95, 25, 82, 10, -31, 32, 46, -6, 61, -3, 93, 40, 61, -9, 35},
        {89, -23, 95, 90, 35, 90, 93, 24, 93, 95, 94, -31, -12, 25, 67, 100, -43, 83, 98, -30, 14, 42, 9, -1, -41, 12, 71, 17, -2, 41, -38, -17},
        {74, 84, 63, 37, 0, 65, 95, -40, 64, 57, 62, 68, 91, 5, -30, 60, 96, -50, 63, 62, 48, 14, -23, 67, -26, -30, -40, 33, 92, 6, 37, 97},
        {72, -11, -23, 8, -17, 48, 58, 51, 91, -42, -49, -18, 37, 60, -26, -3, -27, 97, 33, 78, 24, 43, -37, 91, 76, 86, -33, 32, 22, 34, -8, 38}
    };
    vector<vector<long double>> a(m, vector<long double>(n, 0));
    vector<long double> b(m), c(cost, cost + n);
    for (int i = 0; i < n; i++){
        long double k = ldexpl(i % 2 ? -1 : 1, row_exp[i]);
        for (int j = 0; j < n; j++) a[i][j] = coef[i][j] * k * ldexpl(1, col_exp[j]);
        b[i] = lim[i] * k;
    }
    for (int j = 0; j < n; j++) a[n + j][j] = ldexpl(1, col_exp[j]), b[n + j] = 100;

    BasicSimplex<long double> wide(c, Simplex::MAXIMIZE), dual(b, Simplex::MINIMIZE);
    for (int i = 0; i < m; i++) wide.add_constraint(a[i], b[i], Simplex::LESSEQ);
    for (int j = 0; j < n; j++){
        vector<long double> column(m);
        for (int i = 0; i < m; i++) column[i] = a[i][j];
        dual.add_constraint(column, c[j], Simplex::GREATEQ);
    }
    assert(wide.solve() == Simplex::FEASIBLE);
    assert(dual.solve() == Simplex::FEASIBLE);
    assert(fabsl(wide.value - dual.value) <= 1e-9L * fabsl(wide.value));
    for (int j = 0; j < n; j++){
        long double s = 0, size = fabsl(c[j]);
        for (int i = 0; i < m; i++) s += a[i][j] * dual.val[i], size += fabsl(a[i][j] * dual.val[i]);
        assert(c[j] - s <= 1e-9L * size);
    }
}

/// Nightly seed 124: a 2e-9 roundoff gain in a column with no limiting row made the double dual claim unbounded
void nightly_roundoff_unbounded(){
    int n = 11, m = 13;
    const int c[] = {-10, 9, -18, 24, -9, 16, 6, 1, 18, 11, -14};
    const int ab[13][12] = {
        {3, 5, 0, 3, 5, -4, -2, 6, 1, 2, 0, 16},
        {0, -2, 5, 5, -4, 0, 2, 6, -2, 2, 6, 5},
        {2, 3, 0, 0, -3, 0, 2, -1, 0, -3, -3, 8},
        {0, 3, 4, -2, 0, -3, 2, -2, 5, 0, 0, 18},
        {3, 0, 3, -1, 5, 6, 0, -4, 1, 5, 0, 30},
        {-2, -1, 0, -3, -4, 3, 0, 3, -4, -4, 6, 2},
        {-4, 0, -4, 6, 0, 4, 0, 1, 5, 6, -1, -8},
        {-4, 0, -4, 0, -3, -2, 0, 6, 2, 2, 0, -17},
        {1, 0, -1, 0, -2, 0, 0, -3, -4, 6, 3, -5},
        {2, -4, 2, 0, 0, 0, -1, -2, 0, 5, 2, -6},
        {2, 6, 6, -2, 2, 0, 5, -2, 3, -2, 1, 39},
        {-1, 0, -1, 0, 2, -4, 3, 0, 0, -1, 0, -5},
        {-2, 6, 0, -2, -2, -2, -1, 6, 2, 0, 0, 21}
    };
    vector<double> lim(m);
    for (int i = 0; i < m; i++) lim[i] = ab[i][n];
    Simplex primal(vector<double>(c, c + n), Simplex::MAXIMIZE), dual(lim, Simplex::MINIMIZE);
    for (int i = 0; i < m; i++) primal.add_constraint(vector<double>(ab[i], ab[i] + n), ab[i][n], Simplex::LESSEQ);
    for (int j = 0; j < n; j++){
        vector<double> column(m);
        for (int i = 0; i < m; i++) column[i] = ab[i][j];
        dual.add_constraint(column, c[j], Simplex::GREATEQ);
    }
    assert(primal.solve() == Simplex::FEASIBLE);
    assert(dual.solve() == Simplex::FEASIBLE);
    assert(fabs(primal.value - dual.value) <= 1e-6 * max(1.0, fabs(primal.value)));
}

/// Row scales 2^-13..2^13 and column scales 2^-10..2^10 give variables of 1e5: pruning at EPS broke double past n = 100
template <typename Float>
void scaled_rows_and_columns(long long iterations, int lo, int hi){
    for (long long it = 0; it < iterations; it++){
        int n = stress::rand_int(lo, hi), m = 2 * n;
        vector<long double> c(n), b(m), col(n);
        vector<vector<long double>> a(m, vector<long double>(n, 0));
        for (auto& x : c) x = stress::rand_int(1, 100);
        /// Columns are not rescaled: double broke past column scales of 2^7 apart, long double holds 2^10
        int spread = is_same<Float, double>::value ? 5 : 10;
        for (auto& x : col) x = ldexpl(1, stress::rand_int(-spread, spread));
        for (int i = 0; i < n; i++){
            bool ge = i % 2;
            long double k = ldexpl(ge ? -1 : 1, stress::rand_int(-13, 13));
            for (int j = 0; j < n; j++) a[i][j] = stress::rand_int(-50, 100) * k * col[j];
            b[i] = (ge ? stress::rand_int(1, 500) : stress::rand_int(1000, 10000)) * k;
        }
        for (int j = 0; j < n; j++) a[n + j][j] = col[j], b[n + j] = 100;

        /// Power of two scales keep coefficients exact in both types, so long double and dual optima check the same LP
        BasicSimplex<Float> primal(vector<Float>(c.begin(), c.end()), Simplex::MAXIMIZE);
        BasicSimplex<long double> wide(c, Simplex::MAXIMIZE), dual(b, Simplex::MINIMIZE);
        for (int i = 0; i < m; i++){
            primal.add_constraint(vector<Float>(a[i].begin(), a[i].end()), b[i], Simplex::LESSEQ);
            wide.add_constraint(a[i], b[i], Simplex::LESSEQ);
        }
        for (int j = 0; j < n; j++){
            vector<long double> column(m);
            for (int i = 0; i < m; i++) column[i] = a[i][j];
            dual.add_constraint(column, c[j], Simplex::GREATEQ);
        }
        assert(primal.solve() == Simplex::FEASIBLE);
        assert(wide.solve() == Simplex::FEASIBLE);
        assert(dual.solve() == Simplex::FEASIBLE);

        long double tol = 1e-9L * fabsl(wide.value);
        assert(fabsl(wide.value - dual.value) <= tol);
        /// double stops at reduced costs below its EPS = 1e-9, which left a 1.9e-9 relative gap on seed 101
        long double primal_tol = (is_same<Float, double>::value ? 1e-8L : 1e-9L) * fabsl(wide.value);
        assert(fabsl(primal.value - wide.value) <= primal_tol);
        for (int i = 0; i < m; i++){
            long double s = 0, size = fabsl(b[i]);
            for (int j = 0; j < n; j++) s += a[i][j] * primal.val[j], size += fabsl(a[i][j] * primal.val[j]);
            assert(s - b[i] <= 1e-9L * size);
        }
    }
}

int main(){
    solve_twice_aborts();
    resolve_without_basis_aborts();
    resolve_matches_fresh<double>(stress::scaled(300));
    resolve_matches_fresh<long double>(stress::scaled(100));
    artificial_left_basic();
    uva_1298_rows();
    mixed_magnitude_rows();
    nearly_parallel_rows();
    degenerate_cycling();
    uva_12512_equalities<double>();
    uva_12512_equalities<long double>();
    roundoff_degenerate_rows();
    tiny_ratio_rows(stress::scaled(400));
    nightly_skipped_pivot();
    nightly_roundoff_unbounded();
    scaled_rows_and_columns<double>(stress::scaled(15), 100, 150);
    scaled_rows_and_columns<long double>(stress::scaled(10), 20, 60);
    duality(stress::scaled(300), 1, 40, true);
    duality(stress::scaled(300), 1, 40, false);
    duality(stress::scaled(10), 80, 100, true);
    duality(stress::scaled(30), 50, 100, false);

    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = stress::rand_int(1, it % 5 ? 3 : 4), m = stress::rand_int(0, it % 5 ? 4 : 6);
        Simplex::Sense sense = stress::rand_int(0, 1) ? Simplex::MAXIMIZE : Simplex::MINIMIZE;
        vector<long long> c(n);
        for (auto& x : c) x = stress::rand_int(-5, 5);
        vector<Constraint> cons(m);
        for (auto& k : cons){
            k.a.resize(n);
            for (auto& x : k.a) x = stress::rand_int(0, 2) ? stress::rand_int(-5, 5) : 0;
            k.b = stress::rand_int(-20, 20);
            int t = stress::rand_int(0, 5);
            k.cmp = t < 3 ? Simplex::LESSEQ : t < 5 ? Simplex::GREATEQ : Simplex::EQUAL;
        }

        Fraction small, large;
        bool feasible = brute(n, c, cons, sense, 1000000, small);
        bool feasible_large = brute(n, c, cons, sense, 10000000, large);
        assert(feasible == feasible_large);
        int expected = !feasible ? Simplex::INFEASIBLE : (small.num * large.den != large.num * small.den) ? Simplex::UNBOUNDED : Simplex::FEASIBLE;

        /// Every third LP multiplies each row by 10^k for its own k: the same polytope, but rows of very different size
        Simplex lp(vector<Float>(c.begin(), c.end()), sense);
        for (auto& k : cons){
            long long mult = 1;
            if (it % 3 == 0) for (int e = stress::rand_int(0, 5); e > 0; e--) mult *= 10;
            vector<Float> row(n);
            for (int j = 0; j < n; j++) row[j] = (Float)k.a[j] * mult;
            lp.add_constraint(row, (Float)k.b * mult, k.cmp);
        }

        /// resolve with fresh objectives must agree with the brute force, from the corner the last call left
        auto check_resolves = [&](){
            for (int q = 0; q < 3; q++){
                vector<long long> d(n);
                for (auto& x : d) x = stress::rand_int(-5, 5);
                Fraction lo, hi;
                assert(brute(n, d, cons, sense, 1000000, lo));
                brute(n, d, cons, sense, 10000000, hi);
                int want = lo.num * hi.den != hi.num * lo.den ? Simplex::UNBOUNDED : Simplex::FEASIBLE;
                assert(lp.resolve(vector<Float>(d.begin(), d.end())) == want);
                if (want == Simplex::UNBOUNDED) continue;

                long double got = 0, tol = 1e-6L * max(1.0L, fabsl(lo.value()));
                for (int j = 0; j < n; j++) assert(lp.val[j] >= -1e-7), got += d[j] * lp.val[j];
                assert(fabsl(lp.value - lo.value()) <= tol && fabsl(got - lo.value()) <= tol);
                for (auto& k : cons){
                    long double s = 0;
                    for (int j = 0; j < n; j++) s += k.a[j] * lp.val[j];
                    if (k.cmp != Simplex::GREATEQ) assert(s <= k.b + 1e-6);
                    if (k.cmp != Simplex::LESSEQ) assert(s >= k.b - 1e-6);
                }
            }
        };

        int status = lp.solve();
        Float result = lp.value;
        assert(status == expected);
        if (status == Simplex::UNBOUNDED) check_resolves();
        if (status != Simplex::FEASIBLE) continue;

        long double opt = small.value(), tol = 1e-6L * max(1.0L, fabsl(opt));
        assert(fabsl(result - opt) <= tol);

        long double objective = 0;
        for (int j = 0; j < n; j++){
            assert(lp.val[j] >= -1e-7);
            objective += c[j] * lp.val[j];
        }
        assert(fabsl(objective - opt) <= tol);

        for (auto& k : cons){
            long double s = 0;
            for (int j = 0; j < n; j++) s += k.a[j] * lp.val[j];
            if (k.cmp != Simplex::GREATEQ) assert(s <= k.b + 1e-6);
            if (k.cmp != Simplex::LESSEQ) assert(s >= k.b - 1e-6);
        }
        check_resolves();
    }

    return 0;
}
