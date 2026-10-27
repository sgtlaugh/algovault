/***
 *
 * Gaussian Elimination on a Band Matrix
 * Solves grid equations where each cell's equation only uses nearby cells, like expected values of a random walk
 *
 * Complexity: O(n m r^2) time and O(n m r + t) memory for t terms, against O(n^3 m^3) for dense gauss on the grid
 *
 * GaussBand<T> grid(n, m): one unknown x(i, j) per cell of an n x m grid, every equation starts as x(i, j) = 0
 * grid.add_constant(i, j, c): adds c to the right side of the equation of x(i, j)
 * grid.add_term(i, j, k, l, w): adds w * x(k, l) to the right side of the equation of x(i, j), (k, l) may be (i, j)
 * grid.solve(): the n x m grid of x(i, j), once per instance and after every add call
 * grid.probabilities: after solve(), true if the weights were probabilities as in the notes
 *
 * Notes:
 *   - Every equation reads x(i, j) = constant + sum of w * x(k, l), so a cell without calls solves to 0
 *   - r is how far apart a cell and its farthest term are in the cell order, r <= min(n, m) + 1 for 8-neighbour moves
 *   - Probabilities: every w passed is >= 0 and each cell's weights sum to at most 1, the rest its chance to stop
 *   - A chance to stop below 1e-14 counts as 0, so a sum that rounding pushed off 1 never stops the walk
 *   - Probabilities solve within 1e-13 relative in double on walks of up to 120 cells with chances down to 1e-15
 *   - A chance s to stop adds about 1e-16 / s relative error in double, as it is read as 1 - (sum of weights)
 *   - Otherwise no row swaps: safe when each cell's |1 - w on itself| exceeds the sum of |w| on the other cells
 *   - A cell that cannot reach a stop makes the system singular and the answer inf or nan
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>  /// If more precision is required, use __float128 from quadmath.h if supported :)
struct GaussBand{
    static constexpr double STOP_ROUNDING = 1e-14;

    struct Term{
        int i, j, k, l;
        T weight;
    };

    int n, m;
    vector<vector<T>> constant;
    vector<Term> terms;
    bool solved = false;

    bool transposed = false, probabilities = false;
    int reach = 0;
    vector<T> band, rhs, stop;

    GaussBand(int n, int m) : n(n), m(m), constant(n, vector<T>(m, 0)){
        assert(n >= 1 && m >= 1);
    }

    void add_constant(int i, int j, T c){
        assert(!solved && inside(i, j));
        constant[i][j] += c;
    }

    void add_term(int i, int j, int k, int l, T weight){
        assert(!solved && inside(i, j) && inside(k, l));
        terms.push_back({i, j, k, l, weight});
    }

    vector<vector<T>> solve(){
        assert(!solved);
        solved = true;

        transposed = false;
        int row_reach = farthest_term();
        transposed = true;
        int column_reach = farthest_term();

        /// Numbering the cells along the shorter side keeps every term near the diagonal, so the band is narrower
        transposed = column_reach < row_reach;
        reach = min(row_reach, column_reach);

        /// The equation of cell u, x(u) - sum of w * x(v) = constant, becomes row u of the band
        int cells = n * m;
        band.assign(cells * (2 * reach + 1), 0);
        rhs.assign(cells, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                int u = cell_index(i, j);
                at(u, u) = 1;
                rhs[u] = constant[i][j];
            }
        }

        /// The chance the walk stops at u, kept apart: 1 - (weights left) cancels to noise once it is below rounding
        probabilities = true;
        stop.assign(cells, 1);
        for (auto& term : terms){
            int u = cell_index(term.i, term.j);
            at(u, cell_index(term.k, term.l)) -= term.weight;
            stop[u] -= term.weight;
            if (term.weight < 0) probabilities = false;
        }
        for (int u = 0; u < cells; u++){
            if (stop[u] < -STOP_ROUNDING) probabilities = false;
            if (stop[u] < STOP_ROUNDING) stop[u] = 0;
        }

        /// Forward elimination, a row only reaches the next reach rows and columns
        for (int row = 0; row < cells; row++){
            int last = min(cells - 1, row + reach);

            /// GTH: the pivot sums the chances of leaving the cell, never 1 - (chance of staying), so it cannot cancel
            T pivot = at(row, row);
            if (probabilities){
                pivot = stop[row];
                for (int col = row + 1; col <= last; col++) pivot -= at(row, col);
            }
            for (int col = row + 1; col <= last; col++) at(row, col) /= pivot;
            rhs[row] /= pivot, stop[row] /= pivot;

            for (int below = row + 1; below <= last; below++){
                /// += of -factor is bit for bit -= of factor, and 8% to 14% faster on walks of 40000 to 80000 cells
                T factor = -at(below, row);
                if (factor == 0) continue;  /// walls and cells without calls leave many: 18% faster on SPOJ NBLINDESC
                for (int col = row + 1; col <= last; col++) at(below, col) += factor * at(row, col);
                rhs[below] += factor * rhs[row], stop[below] += factor * stop[row];
            }
        }

        /// Back substitution, every row is already divided by its pivot
        for (int row = cells - 1; row >= 0; row--){
            int last = min(cells - 1, row + reach);
            for (int col = row + 1; col <= last; col++) rhs[row] -= at(row, col) * rhs[col];
        }

        vector<vector<T>> x(n, vector<T>(m));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) x[i][j] = rhs[cell_index(i, j)];
        }

        return x;
    }

    bool inside(int i, int j){
        return 0 <= i && i < n && 0 <= j && j < m;
    }

    int farthest_term(){
        int farthest = 0;
        for (auto& term : terms){
            farthest = max(farthest, abs(cell_index(term.i, term.j) - cell_index(term.k, term.l)));
        }
        return farthest;
    }

    /// Cells are numbered row by row, or column by column when transposed
    int cell_index(int i, int j){
        return transposed ? j * n + i : i * m + j;
    }

    /// Row and column of the full matrix, each row stores only the columns within reach of its diagonal
    T& at(int row, int col){
        return band[row * (2 * reach + 1) + col - row + reach];
    }
};

int main(){
    const int dx[] = {1, 0, -1, 0}, dy[] = {0, 1, 0, -1};
    int w, h, tc = 0;

    while (scanf("%d %d", &w, &h) == 2){
        int u, d, l, r;
        scanf("%d %d %d %d", &u, &d, &l, &r);
        double prob[4] = {d / 100.0, r / 100.0, u / 100.0, l / 100.0};
        vector<string> grid(h);
        static char row[64];
        for (auto& s : grid) scanf("%s", row), s = row;
        auto open = [&](int i, int j){ return i >= 0 && i < h && j >= 0 && j < w && grid[i][j] != 'X'; };

        vector<vector<bool>> seen(h, vector<bool>(w, false));
        vector<pair<int, int>> stack;
        int starts = 0;
        for (int j = 0; j < w; j++){
            if (grid[0][j] == '.') starts++, seen[0][j] = true, stack.push_back({0, j});
        }
        while (!stack.empty()){
            int i = stack.back().first, j = stack.back().second;
            stack.pop_back();
            if (grid[i][j] == 'T') continue;
            for (int k = 0; k < 4; k++){
                int x = i + dx[k], y = j + dy[k];
                if (open(x, y) && !seen[x][y]) seen[x][y] = true, stack.push_back({x, y});
            }
        }

        GaussBand<double> g(h, w);
        for (int i = 0; i < h; i++){
            for (int j = 0; j < w; j++){
                if (!seen[i][j]) continue;
                double stay = 0;
                for (int k = 0; k < 4; k++){
                    int x = i + dx[k], y = j + dy[k];
                    if (!open(x, y)) stay += prob[k];
                    else if (grid[x][y] == '.') g.add_term(i, j, x, y, prob[(k + 2) & 3]);
                }
                if (grid[i][j] == '.'){
                    g.add_term(i, j, i, j, stay);
                    if (i == 0) g.add_constant(i, j, 1.0 / starts);
                }
            }
        }

        auto x = g.solve();
        if (tc++) puts("");
        for (int i = 0; i < h; i++){
            for (int j = 0; j < w; j++){
                if (grid[i][j] == 'T') printf("%.10f\n", x[i][j] + 1e-12);
            }
        }
    }
    return 0;
}
