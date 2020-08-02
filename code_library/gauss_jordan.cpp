/***
Gauss-Jordan Elimination

n = number of linear equations
m = number of variables
equations[i][m] = right-hand side value of constants

For instance, the system of linear equations becomes:

2x + y -z = 8      ----->  (i)
-3x -y + 2z = -11  ----->  (ii)
-2x + y + 2z = -3  ----->  (iii)

n = 3 (x, y, z), m = 3 (i, ii, iii)
equations[0] = {2, 1, -1, 8}    ----->  (i)
equations[1] = {-3, -1, 2, -11} ----->  (ii)
equations[2] = {-2, 1, 2, -3}   ----->  (iii)

Returns -1 when there is no solution
Otherwise returns the number of independent variables, 0 for an unique solution
Contains a valid solution in the vector res on successful completion

Complexity: O(min(n, m) * n * m)

Notes:
For problems on graphs, make sure the graph is connected and a single component
If not, then re-number the vertices and solve for each component separately.

If more precision is required, use long double or __float128 from quadmath.h if supported

***/

#include <bits/stdc++.h>

using namespace std;

template <class T>
int gauss(vector<vector <T>> equations, vector<T>& res, const T eps=1e-12){
    int n = equations.size(), m = equations[0].size() - 1;
    int i, j, k, l, p, f_var = 0;

    res.assign(m, 0);
    vector <int> pos(m, -1);

    for (j = 0, i = 0; j < m && i < n; j++){
        for (k = i, p = i; k < n; k++){
            if (abs(equations[k][j]) > abs(equations[p][j])) p = k;
        }

        if (abs(equations[p][j]) > eps){
            pos[j] = i;
            for (l = j; l <= m; l++) swap(equations[p][l], equations[i][l]);

            for (k = 0; k < n; k++){
                if (k != i){
                    T x = equations[k][j] / equations[i][j];
                    for (l = j; l <= m; l++) equations[k][l] -= equations[i][l] * x;
                }
            }
            i++;
        }
    }

    for (i = 0; i < m; i++){
        if (pos[i] == -1) f_var++;
        else res[i] = equations[pos[i]][m] / equations[pos[i]][i];
    }

    for (i = 0; i < n; i++) {
        T val = 0;
        for (j = 0; j < m; j++) val += res[j] * equations[i][j];
        if (abs(val - equations[i][m]) > eps) return -1;
    }

    return f_var;
}

int main(){
    vector<vector<double>> equations;
    equations = {{2, 1, -1, 8}, {-3, -1, 2, -11}, {-2, 1, 2, -3}};

    vector <double> res;
    int f_var = gauss(equations, res);
    assert(f_var == 0 && res.size() == 3); /// unique solution

    /// x = 2, y = 3, z = -1
    assert(abs(res[0] - 2) < 1e-9);
    assert(abs(res[1] - 3) < 1e-9);
    assert(abs(res[2] + 1) < 1e-9);

    equations = {{2, 1, -1, 8}, {-3, -1, 2, -11}};
    f_var = gauss(equations, res);
    assert(f_var == 1 && res.size() == 3); /// multiple solution, 1 independent variable

    for (int i = 0; i < equations.size(); i++){
        double v = 0;
        int m = res.size();
        for (int j = 0; j < m; j++){
            v += res[j] * equations[i][j];
        }
        assert(abs(v - equations[i][m]) < 1e-9);
    }

    equations = {{2, 1, -1, 8}, {-3, -1, 2, -11}, {4, 2, -2, 10}};
    f_var = gauss(equations, res);
    assert(f_var == -1); /// no solution exists

    return 0;
}
