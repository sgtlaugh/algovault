#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/simplex.cpp"
#undef main

#include <sys/wait.h>
#include <unistd.h>

typedef __int128 i128;
using Float = Simplex::Float;

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

/// Phase 1 enters x0 and row 2 wins the ratio test by 2e-7, leaving the artificial basic at 2e-10, inside double's EPS
/// Left basic it grows in phase 2 and relaxes x0 >= 1000 away, so min x0 came out 0 instead of 1000 or INFEASIBLE
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

    /// Found by random search, each cycles under a different broken tie rule, rows are {coefficients..., rhs} plus sum(x) <= 10
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

/// Primal max c.x, Ax <= b and dual min b.y, A^T y >= c solved as two live instances, strong duality pins their statuses and optima
/// Planted LPs build b and c from x0, y0 >= 0 so both sides are feasible and bounded, arbitrary ones mix in infeasible and unbounded pairs
/// Past n, m = 50 a phase 1 without a ratio test wandered for minutes, the per-instance alarm turns such a hang into a failure
void duality(long long iterations, int lo, int hi, bool planted){
    for (long long it = 0; it < iterations; it++){
        int n = stress::rand_int(lo, hi), m = stress::rand_int(lo, hi);
        vector<vector<long long>> a(m, vector<long long>(n));
        vector<long long> b(m), c(n), x0(n, 0), y0(m, 0);
        for (auto& r : a) for (auto& x : r) x = stress::rand_int(0, 3) ? stress::rand_int(-4, 6) : 0;
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

        alarm(10);
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

/// UVa 1298: three nearly parallel rows of 1e5-sized coefficients; the last row, -6x0 - 3x1 >= 1, cannot hold for x >= 0,
/// yet a ratio test taking the lowest id among near ties pivoted on 5.5e-8 and answered FEASIBLE at x = (0, 0, 1.8e7)
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

/// Rows of very different size: unscaled, the optimum value came out right at x = 3 but y = 1.833173 broke 600000y >= 1100000 by 96
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

/// Nearly parallel rows, infeasible by an exact rational simplex and still infeasible when every limit moves by 1e-7 of its row:
/// rows divided by their largest coefficient but a ratio test taking the lowest id among near ties answered FEASIBLE
void nearly_parallel_rows(){
    Simplex lp({-4, 4, -2}, Simplex::MINIMIZE);
    lp.add_constraint({0, -1, -2}, 3, Simplex::LESSEQ);
    lp.add_constraint({89997, -30003, -90002}, 30000, Simplex::GREATEQ);
    lp.add_constraint({89999, -29997, -90001}, 30000, Simplex::EQUAL);
    lp.add_constraint({90002, -30001, -90002}, 4, Simplex::GREATEQ);
    assert(lp.solve() == Simplex::INFEASIBLE);
}

/// The tableau is pivoted in place, so a second solve() on the same instance must abort instead of answering
void solve_twice_aborts(){
    pid_t pid = fork();
    if (pid == 0){
        close(STDERR_FILENO);
        Simplex encore({3, 2}, Simplex::MAXIMIZE);
        encore.add_constraint({1, 1}, 4, Simplex::LESSEQ);
        encore.solve();
        encore.solve();
        _exit(0);
    }

    int status;
    waitpid(pid, &status, 0);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}

int main(){
    solve_twice_aborts();
    artificial_left_basic();
    uva_1298_rows();
    mixed_magnitude_rows();
    nearly_parallel_rows();
    degenerate_cycling();
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

        int status = lp.solve();
        Float result = lp.value;
        assert(status == expected);
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
    }

    return 0;
}
