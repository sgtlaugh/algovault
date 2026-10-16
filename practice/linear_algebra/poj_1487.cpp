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

/// kind 0 is an integer leaf, 1 an identifier, 2 an inner node
struct Node{
    int kind, val;
    vector<int> ch;
};

vector<Node> nodes;
int n;

int next_char(){
    int c = getchar();
    while (c != EOF && isspace(c)) c = getchar();
    return c;
}

int parse(int c){
    Node node;
    if (c == '('){
        node.kind = 2;
        for (c = next_char(); c != ')'; c = next_char()) node.ch.push_back(parse(c));
    }
    else if (islower(c)){
        node.kind = 1, node.val = c - 'a';
    }
    else{
        node.kind = 0;
        int sign = 1, v = 0;
        if (c == '-') sign = -1, c = getchar();
        for (; isdigit(c); c = getchar()) v = v * 10 + c - '0';
        ungetc(c, stdin);
        node.val = sign * v;
    }
    nodes.push_back(node);
    return nodes.size() - 1;
}

/// Whether some leaf is reachable from node u, given which identifiers reach a leaf
bool reaches_leaf(int u, const vector<bool>& id_leaf){
    if (nodes[u].kind == 0) return true;
    if (nodes[u].kind == 1) return id_leaf[nodes[u].val];
    for (int v : nodes[u].ch){
        if (reaches_leaf(v, id_leaf)) return true;
    }
    return false;
}

/// Every inner node of the tree must reach a leaf, otherwise random play gets stuck below it with positive probability
bool locally_ok(int u, const vector<bool>& id_leaf, vector<bool>& refs){
    if (nodes[u].kind == 1) refs[nodes[u].val] = true;
    if (nodes[u].kind != 2) return true;
    if (!reaches_leaf(u, id_leaf)) return false;
    for (int v : nodes[u].ch){
        if (!locally_ok(v, id_leaf, refs)) return false;
    }
    return true;
}

/// Adds w times the expected score of node u to row, identifiers map to columns through col
void expand(int u, double w, const vector<int>& col, vector<double>& row){
    if (nodes[u].kind == 0) row.back() += w * nodes[u].val;
    else if (nodes[u].kind == 1) row[col[nodes[u].val]] -= w;
    else{
        for (int v : nodes[u].ch) expand(v, w / nodes[u].ch.size(), col, row);
    }
}

int main(){
    for (int game = 1; scanf("%d", &n) == 1 && n; game++){
        nodes.clear();
        vector<int> root(n);
        for (int i = 0; i < n; i++){
            next_char();
            next_char();
            root[i] = parse(next_char());
        }

        vector<bool> id_leaf(n, false);
        for (bool changed = true; changed;){
            changed = false;
            for (int i = 0; i < n; i++){
                if (!id_leaf[i] && reaches_leaf(root[i], id_leaf)) id_leaf[i] = changed = true;
            }
        }

        vector<bool> ok(n);
        vector<vector<bool>> reach(n, vector<bool>(n, false));
        for (int i = 0; i < n; i++){
            ok[i] = locally_ok(root[i], id_leaf, reach[i]);
            reach[i][i] = true;
        }
        for (int k = 0; k < n; k++){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if (reach[i][k] && reach[k][j]) reach[i][j] = true;
                }
            }
        }

        vector<bool> defined(n, true);
        vector<int> col(n, -1), vars;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (reach[i][j] && !ok[j]) defined[i] = false;
            }
            if (defined[i]){
                col[i] = vars.size();
                vars.push_back(i);
            }
        }

        /// A defined identifier only refers to defined ones, so their system alone is nonsingular
        vector<double> res;
        int k = vars.size();
        if (k){
            vector<vector<double>> equations(k, vector<double>(k + 1, 0));
            for (int i = 0; i < k; i++){
                equations[i][i] += 1;
                expand(root[vars[i]], 1, col, equations[i]);
            }
            gauss(equations, res);
        }

        printf("Game %d\n", game);
        for (int i = 0; i < n; i++){
            if (!defined[i]) printf("Expected score for %c undefined\n", 'a' + i);
            else{
                double v = res[col[i]];
                /// Exact answers like 0.9375 land just below the tie after elimination, nudge them away from 0
                if (abs(v) < 1e-9) v = 0;
                else v += v > 0 ? 1e-9 : -1e-9;
                printf("Expected score for %c = %.3f\n", 'a' + i, v);
            }
        }
        puts("");
    }
    return 0;
}
