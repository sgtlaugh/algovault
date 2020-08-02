/***
Gauss-Jordan Elimination modulo primes

n = number of linear equations
m = number of variables
equations[i][m] = right-hand side value of constants

For instance, the system of linear equations becomes:

(2x + y -z) % 7    = 1 % 7  ----->  (i)
(-3x -y + 2z) % 7  = 3 % 7  ----->  (ii)
(-2x + y + 2z) % 7 = 4 % 7  ----->  (iii)

n = 3 (x, y, z), m = 3 (i, ii, iii), mod = 7
equations[0] = {2, 1, 6, 1}  ----->  (i)
equations[1] = {4, 6, 2, 3}  ----->  (ii)
equations[2] = {5, 1, 2, 4}  ----->  (iii)

Returns -1 when there is no solution
Otherwise returns the number of independent variables, 0 for an unique solution
Contains a valid solution in the vector res on successful completion

Complexity: O(min(n, m) * n * m)

Notes:
For problems on graphs, make sure the graph is connected and a single component
If not, then re-number the vertices and solve for each component separately

***/

#include <bits/stdc++.h>

using namespace std;

int expo(int a, int b, int mod){
    int res = 1;

    while (b){
        if (b & 1) res = (long long)res * a % mod;
        a = (long long)a * a % mod;
        b >>= 1;
    }

    return res;
}

/// Note, mod needs to be a prime and equation values should be non-negative
int gauss(vector<vector<int>> equations, vector<int>& res, int mod){
    int n = equations.size(), m = equations[0].size() - 1;

    int i, j, k, l, p, x, inv, f_var = 0;
    const long long modsq = (long long)mod * mod;

    res.assign(m, 0);
    vector <int> pos(m, -1);

    for (j = 0, i = 0; j < m && i < n; j++){
        for (k = i, p = i; k < n; k++){
            if (equations[k][j] > equations[p][j]) p = k;
        }

        if (equations[p][j]){
            pos[j] = i;
            for (l = j; l <= m; l++) swap(equations[p][l], equations[i][l]);

            inv = expo(equations[i][j], mod - 2, mod);
            for (k = 0; k < n && inv; k++){
                if (k != i && equations[k][j]){
                    x = (long long)equations[k][j] * inv % mod;
                    for (l = j; l <= m && x; l++){
                        if (equations[i][l]){
                            equations[k][l] = (modsq + equations[k][l] - (long long)equations[i][l] * x) % mod;
                        }
                    }
                }
            }
            i++;
        }
    }

    for (i = 0; i < m; i++){
        if (pos[i] == -1) f_var++;
        else{
            inv = expo(equations[pos[i]][i], mod - 2, mod);
            res[i] = (long long)equations[pos[i]][m] * inv % mod;
        }
    }

    for (i = 0; i < n; i++) {
        long long val = 0;
        for (j = 0; j < m; j++) val = (val + (long long)res[j] * equations[i][j]) % mod;
        if (val != equations[i][m]) return -1;
    }

    return f_var;
}

int main(){
    const int mod = 7;
    vector<vector<int>> equations;
    equations = {{2, 1, 6, 1}, {4, 6, 2, 3}, {5, 1, 2, 4}};

    vector <int> res;
    int f_var = gauss(equations, res, mod);
    assert(f_var == 0 && res.size() == 3); /// unique solution

    /// x = 2, y = 3, z = -1
    assert(res[0] == 2);
    assert(res[1] == 3);
    assert(res[2] == 6);

    equations = {{2, 1, 6, 1}, {4, 6, 2, 3}};
    f_var = gauss(equations, res, mod);
    assert(f_var == 1 && res.size() == 3); /// multiple solution, 1 independent variable

    for (int i = 0; i < (int)equations.size(); i++){
        long long v = 0;
        int m = res.size();
        for (int j = 0; j < m; j++){
            v = (v + (long long)res[j] * equations[i][j]) % mod;
        }
        assert(v == equations[i][m]);
    }

    equations = {{2, 1, 6, 1}, {4, 6, 2, 3}, {4, 2, 5, 3}};
    f_var = gauss(equations, res, mod);
    assert(f_var == -1); /// no solution exists

    return 0;
}
