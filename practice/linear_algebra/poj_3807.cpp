/***
 *
 * Gauss-Jordan Elimination
 * Linear systems, rank and inverse of matrices over the reals
 *
 * Complexity: O(min(n, m) * n * m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, res): n >= 1 equations in m variables, equations[i][m] is the right-hand side
 *     returns -1 when there is no solution
 *     otherwise returns the number of free variables, 0 for a unique solution, and res holds one solution
 * matrix_rank(a): rank of any n x m matrix
 * matrix_inverse(a, inv): for a square a, returns false when it is singular, otherwise inv = a^-1
 *
 * For instance, the system of linear equations
 *
 *     2x + y - z   = 8    ----->  equations[0] = {2, 1, -1, 8}
 *     -3x - y + 2z = -11  ----->  equations[1] = {-3, -1, 2, -11}
 *     -2x + y + 2z = -3   ----->  equations[2] = {-2, 1, 2, -3}
 *
 * has the unique solution x = 2, y = 3, z = -1
 *
 * For problems on graphs, make sure the graph is connected and a single component
 * If not, then re-number the vertices and solve for each component separately
 *
 * Every function takes an optional eps, scaled by max(1, largest |coefficient| of its input)
 * A pivot within that tolerance of 0 counts as 0, making the matrix rank deficient or singular
 * The right-hand side of gauss is left out of the pivot scale, so x = 1e9 with coefficients of 1 still solves
 * A row of gauss that reduces to 0 = r is consistent when |r| is within eps scaled by the largest |entry| including it
 * If more precision is required, use long double or __float128 from quadmath.h if supported
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Relative, an absolute eps rejects solvable systems with large coefficients
/// Over the first cols columns, a large right-hand side in the scale would push ordinary pivots under it
template <class T>
T tolerance(const vector<vector<T>>& a, int cols, T eps){
    T tol = 1;
    for (auto& row : a) for (int j = 0; j < cols; j++) tol = max(tol, (T)abs(row[j]));
    return tol * eps;
}

/// Eliminates the first cols columns of a with partial pivoting, pivots are not scaled to 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
template <class T>
vector<int> row_reduce(vector<vector<T>>& a, int cols, T tol){
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        for (int k = i; k < n; k++){
            if (abs(a[k][j]) > abs(a[p][j])) p = k;
        }
        if (abs(a[p][j]) <= tol) continue;

        swap(a[p], a[i]);
        int width = a[i].size();
        for (int k = 0; k < n; k++){
            if (k == i || a[k][j] == 0) continue;
            T x = a[k][j] / a[i][j];
            for (int l = j; l < width; l++) a[k][l] -= a[i][l] * x;
        }

        pos[j] = i++;
    }

    return pos;
}

template <class T>
int gauss(vector<vector<T>> equations, vector<T>& res, const T eps=1e-9){
    int n = equations.size(), m = equations[0].size() - 1, f_var = 0;
    T tol = tolerance(equations, m + 1, eps);
    vector<int> pos = row_reduce(equations, m, tolerance(equations, m, eps));

    res.assign(m, 0);
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = nonzero
    for (int k = rank; k < n; k++){
        if (abs(equations[k][m]) > tol) return -1;
    }

    for (int j = 0; j < m; j++){
        if (pos[j] == -1) f_var++;
        else res[j] = equations[pos[j]][m] / equations[pos[j]][j];
    }

    return f_var;
}

template <class T>
bool matrix_inverse(const vector<vector<T>>& a, vector<vector<T>>& inv, const T eps=1e-9){
    int n = a.size();
    vector<vector<T>> aug(n, vector<T>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = row_reduce(aug, n, tolerance(a, n, eps));
    if (n && pos[n - 1] != n - 1) return false;

    inv.assign(n, vector<T>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j] / aug[i][i];
    }

    return true;
}

template <class T>
int matrix_rank(vector<vector<T>> a, const T eps=1e-9){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = row_reduce(a, cols, tolerance(a, cols, eps));
    return cols - count(pos.begin(), pos.end(), -1);
}

typedef long double ld;

map<string, int> element_id;

int read_number(const string& s, size_t& i){
    if (i >= s.size() || !isdigit(s[i])) return 1;
    int v = 0;
    while (i < s.size() && isdigit(s[i])) v = v * 10 + s[i++] - '0';
    return v;
}

int read_element(const string& s, size_t& i){
    string e(1, s[i++]);
    if (i < s.size() && islower(s[i])) e += s[i++];
    if (!element_id.count(e)){
        int k = element_id.size();
        element_id[e] = k;
    }
    return element_id[e];
}

map<int, int> molecule(const string& s){
    map<int, int> cnt;
    for (size_t i = 0; i < s.size();){
        if (s[i] == '('){
            map<int, int> group;
            for (i++; s[i] != ')';){
                int e = read_element(s, i);
                group[e] += read_number(s, i);
            }
            int mul = read_number(s, ++i);
            for (auto& g : group) cnt[g.first] += g.second * mul;
        }
        else{
            int e = read_element(s, i);
            cnt[e] += read_number(s, i);
        }
    }
    return cnt;
}

int main(){
    string line;
    while (cin >> line && line != "."){
        line.pop_back();
        size_t arrow = line.find("->");
        element_id.clear();

        vector<map<int, int>> mols;
        vector<int> sign;
        for (int side = 0; side < 2; side++){
            string part = side ? line.substr(arrow + 2) : line.substr(0, arrow);
            stringstream ss(part);
            string mol;
            while (getline(ss, mol, '+')){
                mols.push_back(molecule(mol));
                sign.push_back(side ? -1 : 1);
            }
        }

        /// The solutions form a line, pinning the last coefficient to 1 makes the solution unique
        int m = mols.size(), e = element_id.size();
        vector<vector<ld>> equations(e + 1, vector<ld>(m + 1, 0));
        for (int j = 0; j < m; j++){
            for (auto& c : mols[j]) equations[c.first][j] = sign[j] * c.second;
        }
        equations[e][m - 1] = equations[e][m] = 1;

        vector<ld> res;
        gauss(equations, res);

        /// The true last coefficient is the smallest scale turning every ratio into an integer
        for (int last = 1;; last++){
            bool integral = true;
            for (int j = 0; j < m && integral; j++){
                ld v = res[j] * last;
                if (abs(v - roundl(v)) > 1e-6) integral = false;
            }
            if (integral){
                for (int j = 0; j < m; j++) printf("%lld%c", (long long)roundl(res[j] * last), j + 1 < m ? ' ' : '\n');
                break;
            }
        }
    }
    return 0;
}
