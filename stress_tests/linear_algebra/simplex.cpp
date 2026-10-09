#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/simplex.cpp"
#undef main

#include <unistd.h>

typedef __int128 i128;

struct Constraint{
    vector<long long> a;
    long long b;
    int cmp;
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
bool brute(int n, const vector<long long>& c, const vector<Constraint>& cons, int sense, long long box, Fraction& best){
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
                if ((k.cmp == LESSEQ && s > rhs) || (k.cmp == GREATEQ && s < rhs) || (k.cmp == EQUAL && s != rhs)) return;
            }

            Fraction value{0, d};
            for (int j = 0; j < n; j++) value.num += c[j] * num[j];
            if (sense == MINIMIZE) value.num = -value.num;
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
    if (found && sense == MINIMIZE) best.num = -best.num;
    return found;
}

/// Degenerate LPs that cycle forever without anti-cycling, the alarm turns a hang into a failure
void degenerate_cycling(){
    Float result;
    alarm(2);

    /// Beale (1955): Dantzig's rule with first-row ratio ties cycles here, optimum 1/20 at x1 = 1/25, x3 = 1
    Simplex::init(4, (Float[]){0, 0.75, -150, 0.02, -6}, MAXIMIZE);
    Simplex::add_constraint((Float[]){0, 0.25, -60, -0.04, 9}, 0, LESSEQ);
    Simplex::add_constraint((Float[]){0, 0.5, -90, -0.02, 3}, 0, LESSEQ);
    Simplex::add_constraint((Float[]){0, 0, 0, 1, 0}, 1, LESSEQ);
    assert(Simplex::solve(result) == FEASIBLE);
    assert(abs(result - 0.05) < 1e-9);
    assert(abs(Simplex::val[1] - 0.04) < 1e-9 && abs(Simplex::val[3] - 1) < 1e-9);

    /// Found by random search, each cycles under a different broken tie rule, rows are {coefficients..., rhs} plus sum(x) <= 10
    vector<tuple<vector<Float>, vector<vector<Float>>, Float>> cases = {
        {{3, -2, 2, 6, 6}, {{-2, 0, 0, -1, -3, 0}, {0, 3, 0, 0, 0, 0}, {0, 3, 4, 0, -3, 1}, {0, -4, 1, 4, 2, 0}, {0, -2, -1, -2, 4, 1}, {1, 0, 0, 1, 3, 0}, {4, 2, 2, 0, 2, 0}}, 0},
        {{3, -1, -1, -1, -4, 0, -2}, {{0, -4, 0, 0, 2, 1, 2, 0}, {4, 3, 0, -1, 0, -1, 1, 0}, {2, 0, -2, -4, 2, 0, 2, 0}, {-2, 0, 0, 4, 3, -2, 0, 0}, {0, 4, 1, 3, 3, -2, 2, 0}, {1, 0, 2, 0, 0, -4, 3, 0}, {-3, -4, -3, -3, 0, 0, 0, 0}}, 0},
        {{-1, 5, -3, -1, -1, 3, 2}, {{-2, 0, -3, 1, 1, 2, 1, 3}, {0, 2, 2, -2, 0, -2, 0, 0}, {3, 4, -2, 4, 0, 4, -3, 0}, {2, 1, -3, 0, 3, 0, -4, 0}, {0, -2, 4, 2, 0, 0, 0, 0}, {-3, 0, 1, 2, 1, 0, 0, 0}, {0, -1, -3, 0, 0, -4, -3, 0}}, 123.0 / 7},
        {{-3, 2, -1, -2, -3, 3, 4}, {{-4, -4, 0, -1, -2, -1, 4, 0}, {-1, 2, 1, 0, -2, 2, -2, 0}, {0, 2, 3, -1, 0, 0, -1, 0}, {4, 0, -2, -3, -3, 0, 3, 0}, {2, -3, 0, 1, -2, 1, 4, 0}}, 65.0 / 8},
    };
    for (auto& [obj, rows, opt] : cases){
        int n = obj.size();
        vector<Float> row(n + 1, 0);
        copy(obj.begin(), obj.end(), row.begin() + 1);
        Simplex::init(n, row.data(), MAXIMIZE);
        for (auto& r : rows){
            copy(r.begin(), r.begin() + n, row.begin() + 1);
            Simplex::add_constraint(row.data(), r[n], LESSEQ);
        }
        fill(row.begin() + 1, row.end(), 1);
        Simplex::add_constraint(row.data(), 10, LESSEQ);
        assert(Simplex::solve(result) == FEASIBLE);
        assert(abs(result - opt) < 1e-9);
    }
    alarm(0);
}

int main(){
    degenerate_cycling();
    srand(stress::seed());
    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = stress::rand_int(1, it % 5 ? 3 : 4), m = stress::rand_int(0, it % 5 ? 4 : 6), sense = stress::rand_int(0, 1) ? MAXIMIZE : MINIMIZE;
        vector<long long> c(n);
        for (auto& x : c) x = stress::rand_int(-5, 5);
        vector<Constraint> cons(m);
        for (auto& k : cons){
            k.a.resize(n);
            for (auto& x : k.a) x = stress::rand_int(0, 2) ? stress::rand_int(-5, 5) : 0;
            k.b = stress::rand_int(-20, 20);
            int t = stress::rand_int(0, 5);
            k.cmp = t < 3 ? LESSEQ : t < 5 ? GREATEQ : EQUAL;
        }

        Fraction small, large;
        bool feasible = brute(n, c, cons, sense, 1000000, small);
        bool feasible_large = brute(n, c, cons, sense, 10000000, large);
        assert(feasible == feasible_large);
        int expected = !feasible ? INFEASIBLE : (small.num * large.den != large.num * small.den) ? UNBOUNDED : FEASIBLE;

        vector<Float> obj(n + 1, 0);
        for (int j = 0; j < n; j++) obj[j + 1] = c[j];
        Simplex::init(n, obj.data(), sense);
        for (auto& k : cons){
            vector<Float> coef(n + 1, 0);
            for (int j = 0; j < n; j++) coef[j + 1] = k.a[j];
            Simplex::add_constraint(coef.data(), k.b, k.cmp);
        }

        Float result;
        int status = Simplex::solve(result);
        assert(status == expected);
        if (status != FEASIBLE) continue;

        long double opt = small.value(), tol = 1e-6L * max(1.0L, fabsl(opt));
        assert(fabsl(result - opt) <= tol);

        long double objective = 0;
        for (int j = 0; j < n; j++){
            assert(Simplex::val[j + 1] >= -1e-7);
            objective += c[j] * Simplex::val[j + 1];
        }
        assert(fabsl(objective - opt) <= tol);
        for (auto& k : cons){
            long double s = 0;
            for (int j = 0; j < n; j++) s += k.a[j] * Simplex::val[j + 1];
            if (k.cmp != GREATEQ) assert(s <= k.b + 1e-6);
            if (k.cmp != LESSEQ) assert(s >= k.b - 1e-6);
        }
    }
    return 0;
}
