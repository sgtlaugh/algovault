/***
 *
 * Simplex Algorithm for Linear Programming
 * Solves maximize/minimize linear objective subject to linear constraints
 *
 * Complexity: O(2^n) worst case, polynomial average case
 * Optimized with sparse pivot - only updates non-zero coefficients
 *
 * Usage:
 *   1. init(n, objective[], MAXIMIZE or MINIMIZE)
 *   2. add_constraint(coefficients[], rhs, LESSEQ or GREATEQ or EQUAL)
 *   3. solve(result) returns FEASIBLE, INFEASIBLE, or UNBOUNDED
 *   4. If feasible, result has optimal value, val[] has variable values
 *
 * Variables x1, x2, ..., xn (1-indexed) satisfy xi >= 0
 * For unbounded variable: replace x with (x1 - x2) where x1, x2 >= 0
 * For constraint |x| <= M: add both x <= M and -x <= M
 *
 * Change `using Float = long double;` for higher precision (auto-adjusts EPS)
 *
 * Example:
 *   Maximize 3x1 + 2x2 subject to x1 + x2 <= 4, 2x1 + x2 <= 5
 *
 *   Float result;
 *   Simplex::init(2, (Float[]){0, 3, 2}, MAXIMIZE);
 *   Simplex::add_constraint((Float[]){0, 1, 1}, 4, LESSEQ);
 *   Simplex::add_constraint((Float[]){0, 2, 1}, 5, LESSEQ);
 *   int status = Simplex::solve(result);
 * 
 *   // result = 9, val[1] = 1, val[2] = 3
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

using Float = double;  // Change to 'long double' for higher precision

const Float EPS = std::is_same<Float, long double>::value ? 1e-13L : 1e-9;

const int MAXV = 4010;
const int MAXC = 4010;

const int EQUAL = 0;
const int LESSEQ = -1;
const int GREATEQ = 1;
const int MINIMIZE = -1;
const int MAXIMIZE = 1;

const int FEASIBLE = 1;
const int INFEASIBLE = -1;
const int UNBOUNDED = 666;

namespace Simplex {
    int n, m, flag;
    Float ar[MAXC][MAXV], val[MAXV], rhs[MAXC];
    int adj[MAXV], idx[MAXV], down[MAXV], link[MAXC];

    void init(int nvar, Float obj[], int min_or_max){
        m = 0, n = nvar, flag = min_or_max;
        memset(ar, 0, sizeof(ar));
        memset(val, 0, sizeof(val));
        memset(rhs, 0, sizeof(rhs));
        memset(idx, 0, sizeof(idx));
        for (int i = 1; i <= n; i++) ar[0][i] = obj[i] * flag;
    }

    void add_constraint(Float coef[], Float lim, int cmp){
        cmp *= -1;
        if (cmp == 0){
            rhs[++m] = lim;
            for (int i = 1; i <= n; i++) ar[m][i] = coef[i];
            rhs[++m] = -lim;
            for (int i = 1; i <= n; i++) ar[m][i] = -coef[i];
        }
        else{
            rhs[++m] = lim * cmp;
            for (int i = 1; i <= n; i++) ar[m][i] = coef[i] * cmp;
        }
    }

    void pivot(int x, int y, Float& res){
        int i, j, len = 0;
        Float v = ar[x][y];

        swap(link[x], down[y]);
        rhs[x] /= v, ar[x][y] = 1;
        for (j = 1; j <= n; j++){
            ar[x][j] /= v;
            if (abs(ar[x][j]) > EPS) adj[len++] = j;
        }

        for (i = 1; i <= m; i++){
            if (abs(ar[i][y]) > EPS && i != x){
                rhs[i] -= ar[i][y] * rhs[x], v = ar[i][y], ar[i][y] = 0;
                for (j = 0; j < len; j++) ar[i][adj[j]] -= (v * ar[x][adj[j]]);
            }
        }

        res += (ar[0][y] * rhs[x]), v = ar[0][y], ar[0][y] = 0;
        for (j = 0; j < len; j++) ar[0][adj[j]] -= (v * ar[x][adj[j]]);
    }

    int solve(Float& res){
        res = 0;
        int i, j, x, y;
        Float u, v, mn, mx;

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
        while (true){
            x = 0, y = 0, mx = EPS;
            for (i = 1; i <= n; i++){
                if (ar[0][i] > mx) mx = ar[0][i], y = i;
            }
            if (y == 0) break;

            for (i = 1; i <= m; i++){
                if (ar[i][y] > EPS){
                    u = rhs[i] / ar[i][y];
                    if (x == 0 || u < v) x = i, v = u;
                }
            }
            if (x == 0) return UNBOUNDED;
            pivot(x, y, res);
        }

        res *= flag;
        for (int i = 1; i <= m; i++){
            if (link[i] <= n) idx[link[i]] = i;
        }

        for (int i = 1; i <= n; i++) val[i] = rhs[idx[i]];
        return FEASIBLE;
    }
}

int main(){
	// Production planning problem
    // Maximize profit: 5x1 + 4x2 + 3x3
    // Subject to:      2x1 + 3x2 + x3 <= 100  (material A)
    //                  4x1 + x2 + 2x3 <= 80   (material B)
    //                  3x1 + 4x2 + 2x3 <= 120 (labor hours)
    
    Float result = 0;

    Simplex::init(3, (Float[]){0, 5, 4, 3}, MAXIMIZE);
    Simplex::add_constraint((Float[]){0, 2, 3, 1}, 100, LESSEQ);
    Simplex::add_constraint((Float[]){0, 4, 1, 2}, 80, LESSEQ);
    Simplex::add_constraint((Float[]){0, 3, 4, 2}, 120, LESSEQ);

    assert(Simplex::solve(result) == FEASIBLE);
    assert(abs(result - 153.333333333) < 1e-6);

    // Infeasible constraints
    Simplex::init(2, (Float[]){0, 1, 1}, MAXIMIZE);
    Simplex::add_constraint((Float[]){0, 1, 1}, 10, LESSEQ);  // x1 + x2 <= 10
    Simplex::add_constraint((Float[]){0, 1, 1}, 20, GREATEQ); // x1 + x2 >= 20 (contradiction!)

    assert(Simplex::solve(result) == INFEASIBLE);

    return 0;
}
