/***
 *
 * Thomas Algorithm
 * Solves a tridiagonal system of n linear equations by forward elimination and back substitution
 *
 * Complexity: O(n) time and memory
 *
 * equation(l, p, r, rhs): the row (x_prev * l) + (x_cur * p) + (x_next * r) = rhs, l of the first row and r of the last are ignored
 * thomas_algorithm(ar): returns x_0..x_{n-1} for the n = ar.size() rows, no pivoting so it needs a diagonally dominant system (|p| >= |l| + |r|), a zero pivot gives inf or nan
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct equation{
    long double l, p, r, rhs;

    equation() {}
    equation(long double l, long double p, long double r, long double rhs = 0.0):
        l(l), p(p), r(r), rhs(rhs) {}
};

vector<long double> thomas_algorithm(vector<struct equation> ar){
    int n = ar.size();
    ar[0].r = ar[0].r / ar[0].p;
    ar[0].rhs = ar[0].rhs / ar[0].p;

    for (int i = 1; i < n; i++){
        long double v = 1.0 / (ar[i].p - ar[i].l * ar[i - 1].r);
        ar[i].r = ar[i].r * v;
        ar[i].rhs = (ar[i].rhs - ar[i].l * ar[i - 1].rhs) * v;
    }

    for (int i = n - 2; i >= 0; i--) ar[i].rhs = ar[i].rhs - ar[i].r * ar[i + 1].rhs;

    vector<long double> res;
    for (int i = 0; i < n; i++) res.push_back(ar[i].rhs);
    return res;
}

int main(){
    /// 2x + y = 4, x + 3y + z = 9, y + 2z = 6 has the solution x = 1, y = 2, z = 2
    auto res = thomas_algorithm({equation(0, 2, 1, 4), equation(1, 3, 1, 9), equation(1, 2, 0, 6)});
    assert(res.size() == 3 && fabsl(res[0] - 1) < 1e-12 && fabsl(res[1] - 2) < 1e-12 && fabsl(res[2] - 2) < 1e-12);

    /// 4x = 2, a lone row ignores both l and r
    res = thomas_algorithm({equation(7, 4, 9, 2)});
    assert(res.size() == 1 && fabsl(res[0] - 0.5) < 1e-12);

    return 0;
}
