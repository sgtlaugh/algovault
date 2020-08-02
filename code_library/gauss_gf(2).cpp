/***
Gauss-Jordan Elimination over GF(2)
Can be done with standard gauss but we can use bitsets to make it faster

n = number of linear equations
m = number of variables
ar[i][m] = right-hand side value of constants

For instance, the system of linear equations (note not in GF(2)) becomes:

(2x + y -z) % 2    = 0 % 2  ----->  (i)
(-3x -y + 2z) % 2  = 1 % 2  ----->  (ii)
(-2x + y + 2z) % 2 = 1 % 2  ----->  (iii)

n = 3 (x, y, z), m = 3 (i, ii, iii)
ar[0] = {0, 1, 1, 0}  ----->  (i)
ar[1] = {1, 1, 0, 1}  ----->  (ii)
ar[2] = {0, 1, 0, 1}  ----->  (iii)

Complexity: O(min(n, m) * n * m / 32)
It's actually O(min(n, m) * n * MAX / 32) because of the implementation since C++ bitset needs to be of constant size

Returns -1 when there is no solution
Otherwise returns the number of independent variables, 0 for an unique solution
Contains a solution in the bit vector res on successful completion

To make it run even faster in practice, use custom bitset class with unsigned long long

***/

#include <bits/stdc++.h>

using namespace std;

const int MAX = 1024;

int gauss(int m, vector <bitset<MAX>> equations, bitset <MAX>& res){
    res.reset();
    vector <int> pos(m, -1);
    int i, j, k, l, v, p, n = equations.size(), f_var = 0;

    for (j = 0, i = 0; j < m && i < n; j++){
        for (p = i, k = i + 1; p == i && k < n; k++){
            if (equations[k][j]) p = k;
        }

        if (equations[p][j]){
            pos[j] = i;
            swap(equations[p], equations[i]);

            for (k = 0; k < n; k++){
                if (k != i && equations[k][j]) equations[k] ^= equations[i];
            }
            i++;
        }
    }

    for (i = 0; i < m; i++){
        if (pos[i] == -1) f_var++;
        else res[i] = equations[pos[i]][m];
    }

    for (i = 0; i < n; i++) {
        for (j = 0, v = 0; j < m; j++) v ^= (res[j] & equations[i][j]);
        if (v != equations[i][m]) return -1;
    }

    return f_var;
}

int main(){
    bitset <MAX> equation;
    vector <bitset<MAX>> equations;

    equation.reset();
    equation[0] = 0, equation[1] = 1, equation[2] = 1, equation[3] = 0;
    equations.push_back(equation);

    equation.reset();
    equation[0] = 1, equation[1] = 1, equation[2] = 0, equation[3] = 1;
    equations.push_back(equation);

    equation.reset();
    equation[0] = 0, equation[1] = 1, equation[2] = 0, equation[3] = 1;
    equations.push_back(equation);

    const int m = 3; /// 3 variables, 1 RHS constant

    bitset <MAX> res;
    int f_var = gauss(m, equations, res);
    assert(f_var == 0);

    assert(res[0] == 0);
    assert(res[1] == 1);
    assert(res[2] == 1);

    return 0;
}
