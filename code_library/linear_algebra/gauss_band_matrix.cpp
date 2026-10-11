/***
 *
 * Gaussian Elimination on a Band Matrix
 * Solves the linear system of a grid where each cell's equation only touches nearby cells
 *
 * Complexity: O(n m band_size^2) time and O(n m band_size) memory
 * With the default band_size = 2 m + 3 that is O(n m^3) time and O(n m^2) memory, vs O(n^3 m^3) for standard gauss on the grid
 *
 * GaussBand<T> g(n, m, rhs = 0, band_size = 0): one unknown x(i, j) and one equation per cell of the n x m grid
 *   every equation starts as 1 * x(i, j) = rhs
 * g.add_coef(i, j, k, l, c): adds c to the coefficient of x(k, l) in the equation of (i, j), so repeated moves accumulate
 *   (k, l) must lie within band_size / 2 of (i, j) in row-major order, asserted; the default reaches (i +- 1, j +- 1)
 * g.set_rhs(i, j, val): replaces the right-hand side of the equation of (i, j)
 * g.solve(): the n x m grid of x(i, j), eliminates in place so it runs once per instance (asserted)
 *
 * Pivots on the diagonal without row swaps, so it needs every pivot to stay nonzero
 * Diagonally dominant systems (such as expected value equations) are safe, others can give inf or nan
 * If n < m, transpose the grid: the band grows with m
 *
 * Expected moves of a random walk to the corner (n - 1, m - 1), each move goes to one of the d neighbours:
 *   GaussBand<double> g(n, m, 1);                        // e(i, j) - sum e(next) / d = 1
 *   for each cell (i, j) but the corner, for each neighbour (k, l): g.add_coef(i, j, k, l, -1.0 / d);
 *   g.set_rhs(n - 1, m - 1, 0);                          // e(corner) = 0
 *   auto e = g.solve();
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>  /// If more precision is required, use __float128 from quadmath.h if supported :)
struct GaussBand{
    int n, m, h, band_size;

    GaussBand(int n, int m, T rhs_default=0, int band_size=0) : n(n), m(m), band_size(band_size ? band_size : 2 * m + 3){
        assert(n >= 1 && m >= 1 && this->band_size >= 1 && this->band_size % 2 == 1);
        h = this->band_size / 2;

        band.assign(n * m * this->band_size, 0);
        for (int i = 0; i < n * m; i++) band[i * this->band_size + h] = 1;
        rhs.assign(n * m, rhs_default);
    }

    void add_coef(int i, int j, int k, int l, T val){
        assert(0 <= i && i < n && 0 <= j && j < m && 0 <= k && k < n && 0 <= l && l < m);
        int u = get_cell_num(i, j), v = get_cell_num(k, l) - u + h;
        assert(0 <= v && v < band_size);
        band[u * band_size + v] += val;
    }

    void set_rhs(int i, int j, T val){
        rhs[get_cell_num(i, j)] = val;
    }

    vector<vector<T>> solve(){
        assert(!band.empty());
        gauss();
        band.clear();  /// the eliminated band is useless and an empty band makes a second solve() assert

        vector<vector<T>> res(n, vector<T>(m));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) res[i][j] = rhs[get_cell_num(i, j)];
        }

        return res;
    }

private:
    vector<T> band, rhs;

    /// Cells are numbered in reverse row-major order, the band offsets of a cell's neighbours are what matter
    int get_cell_num(int i, int j){
        return m * (n - i) - j - 1;
    }

    void gauss(){
        T x;
        int i, j, k, l, d, u, v, q = n * m;

        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                d = i * m + j;
                x = band[d * band_size + h], rhs[d] /= x;
                for (k = 0; k <= h && (d + k) < q; k++) band[d * band_size + h + k] /= x;

                for (l = 1; l <= h && (d + l) < q; l++){
                    x = band[(d + l) * band_size + h - l], rhs[d + l] -= x * rhs[d];
                    for (k = 0; k <= h && (d + k) < q; k++){
                        band[(d + l) * band_size + h - l + k] -= x * band[d * band_size + h + k];
                    }
                }
            }
        }

        vector<T> aux(q, 0);
        swap(aux, rhs);
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                u = get_cell_num(i, j);
                for (x = 0, v = 0; v < band_size; v++){
                    int idx = u + v - h;
                    if (idx >= 0 && idx < n * m) x += band[u * band_size + v] * rhs[idx];
                }
                rhs[u] = aux[u] - x;
            }
        }
    }
};

int main(){
    /// Random walk on a 2 x 2 grid to (1, 1): e(0, 0) = 1 + (e(0, 1) + e(1, 0)) / 2, e(0, 1) = e(1, 0) = 1 + (e(0, 0) + e(1, 1)) / 2
    GaussBand<long double> walk(2, 2, 1);
    walk.set_rhs(1, 1, 0);
    walk.add_coef(1, 0, 0, 0, -0.5), walk.add_coef(1, 0, 1, 1, -0.5);
    walk.add_coef(0, 1, 0, 0, -0.5), walk.add_coef(0, 1, 1, 1, -0.5);
    walk.add_coef(0, 0, 1, 0, -0.5), walk.add_coef(0, 0, 0, 1, -0.5);

    auto e = walk.solve();
    assert(abs(e[0][0] - 4) < 1e-9 && abs(e[0][1] - 3) < 1e-9 && abs(e[1][0] - 3) < 1e-9 && abs(e[1][1]) < 1e-9);

    /// 1 x 3 corridor to (0, 2), a move off either end stays put: e(0, 0) = 1 + (e(0, 0) + e(0, 1)) / 2 gives e = {6, 4, 0}
    GaussBand<double> corridor(1, 3, 1);
    corridor.set_rhs(0, 2, 0);
    corridor.add_coef(0, 0, 0, 0, -0.5), corridor.add_coef(0, 0, 0, 1, -0.5);
    corridor.add_coef(0, 1, 0, 0, -0.5), corridor.add_coef(0, 1, 0, 2, -0.5);

    auto c = corridor.solve();
    assert(abs(c[0][0] - 6) < 1e-9 && abs(c[0][1] - 4) < 1e-9 && abs(c[0][2]) < 1e-9);

    return 0;
}
