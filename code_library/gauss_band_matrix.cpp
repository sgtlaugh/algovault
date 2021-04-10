/***
 * Gaussian Elimination on sparse band matrix
 * Used to solve systems of linear equations in grids
 * Usually equations are formed by linear combinations of variables in the same/adjcent rows
 * Since not all variables are used we can group them in "bands" and solve faster
 *
 * Standard gauss in grids complexity: O(n^3 * m^3)
 * Gauss with band matrix: O(n * m^3) or O(n * m * band_size ^ 2) to be more precise
 * If n < m, better to transpose the matrix
 *
***/

#include <bits/stdc++.h>

#define MAX_ROW  10010
#define MAX_COL  32
#define MAX_BAND 70

using namespace std;

template <class T>  /// If more precision is required, use __float128 from quadmath.h if supported :)
struct GaussBand{
    int n, m, h, band_size;
    T band[MAX_ROW * MAX_COL][MAX_BAND];
    T rhs[MAX_ROW * MAX_COL], aux[MAX_ROW * MAX_COL];

    /***
     * Initialize a gauss band matrix
     * n = number of rows
     * m = number of columns
     * rhs_default is the default RHS constant for all equations, can be explicitly updated with set_rhs() for specific cells
     * band_size = maximum band size, the equation for any variable is a linear combination of only variables inside the band
     *
    ***/
    GaussBand(int n, int m, int band_size=0, T rhs_default=0) : n(n), m(m) {
        if (!band_size) band_size = 2 * m + 3;
        assert(n >= 1 && n < MAX_ROW && m >= 1 && m < MAX_COL
               && band_size >= 1 && band_size < MAX_BAND && band_size % 2 == 1);

        this->band_size = band_size;
        this->h = band_size / 2;

        for (int i = 0; i < n * m; i++){
            for (int j = 0; j < band_size; j++){
                band[i][j] = 0;
            }
            band[i][h] = 1;
            rhs[i] = rhs_default;
        }
    }

    /***
     * Transforms 2D array to 1D for easier intermediate calculation
     *
    ***/
    int get_cell_num(int i, int j){
        return m * (n - i) - j - 1;
    }

    /***
     * Sets gauss_matrix[u][v] to val, where u = get_cell_num(i, j) and v = get_cell_num(k, l)
     * (i, j) and (k, l) are cells in the original matrix
     * gauss_matrix[u] contains the coefficients of the equation for cell(i, j)
     * gauss_matrix does not exist, the coordinates are converted to band matrix implicitly
     *
    ***/
    void set_matrix(int i, int j, int k, int l, T val){
        int u = get_cell_num(i, j), v = get_cell_num(k, l) - u + h;
        band[u][v] = val;
    }

    /***
     * Sets rhs[u] to val, where u = get_cell_num(i, j)
     * rhs[u] contains the right hand constant for equation[u]
     *
    ***/
    void set_rhs(int i, int j, T val){
        rhs[get_cell_num(i, j)] = val;
    }

    void gauss(){
        T x;
        int i, j, k, l, d, u, v, q = n * m;

        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                d = i * m + j;
                x = band[d][h], rhs[d] /= x;
                for (k = 0; k <= h && (d + k) < q; k++) band[d][h + k] /= x;

                for (l = 1; l <= h && (d + l) < q; l++){
                    x = band[d + l][h - l], rhs[d + l] -= x * rhs[d];
                    for (k = 0; k <= h && (d + k) < q; k++){
                        band[d + l][h - l + k] -= x * band[d][h + k];
                    }
                }
            }
        }

        for (i = 0; i < q; i++) aux[i] = rhs[i], rhs[i] = 0;
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                u = get_cell_num(i, j);
                for (x = 0, v = 0; v < band_size; v++){
                    int idx = u + v - h;
                    if (idx >= 0 && idx < n * m) x += band[u][v] * rhs[idx];
                }
                rhs[u] = aux[u] - x;
            }
        }
    }

    vector <vector<T>> solve(){
        gauss();

        vector <vector<T>> res;
        for (int i = 0; i < n; i++){
            vector <T> v;
            for (int j = 0; j < m; j++){
                int idx = get_cell_num(i, j);
                v.push_back(rhs[idx]);
            }
            res.push_back(v);
        }

        return res;
    }
};

int main(){
    /***
     * Let's solve the following problem
     * In a n x m grid, player starts from a starting cell
     * From each cell, moves randomly to any adjacent valid cell with equal probability
     * Journey ends when cell(n - 1, m - 1) is reached
     * Find the expected number of moves for each starting cell
     *
     * Let's write the equations for n=2, m=2,
     * e(i, j) is the expected number of moves to end from cell(i, j)
     *
     * e(1, 1) = 0
     * e(1, 0) = 0.5 * e(0, 0) + 0.5 * e(1, 1) + 1
     * e(0, 1) = 0.5 * e(0, 0) + 0.5 * e(1, 1) + 1
     * e(0, 0) = 0.5 * e(0, 1) + 0.5 * e(1, 0) + 1
     *
     * Or, after re-arranging:
     * e(1, 1) = 0
     * e(1, 0) - 0.5 * e(0, 0) - 0.5 * e(1, 1) = 1
     * e(0, 1) - 0.5 * e(0, 0) - 0.5 * e(1, 1) = 1
     * e(0, 0) - 0.5 * e(0, 1) - 0.5 * e(1, 0) = 1
     *
    ***/

    int n = 2, m = 2;

    /***
     * Write all the cells in row major order and number them
     * For any cell(i, j), the maximum band size is calculated as below
     * The difference between the numbers in immediate top and bottom cell
     * Or |num(i - 1, j) - num(i + 1, j)| + 1, in this case 2*m + 1
     *
    ***/
    int band_size = 2 * m + 1;

    /// RHS constant for most variables except (1, 1) is 1
    long double rhs_default = 1;

    auto gauss_band = GaussBand <long double>(n, m, band_size=band_size, rhs_default=rhs_default);

    gauss_band.set_rhs(1, 1, 0);

    gauss_band.set_matrix(1, 0, 0, 0, -0.5);  /// coefficient of e(0, 0) in the equation for e(1, 0)
    gauss_band.set_matrix(1, 0, 1, 1, -0.5);  /// coefficient of e(1, 1) in the equation for e(1, 0)

    gauss_band.set_matrix(0, 1, 0, 0, -0.5);  /// coefficient of e(0, 0) in the equation for e(0, 1)
    gauss_band.set_matrix(0, 1, 1, 1, -0.5);  /// coefficient of e(1, 1) in the equation for e(0, 1)

    gauss_band.set_matrix(0, 0, 1, 0, -0.5);  /// coefficient of e(1, 0) in the equation for e(0, 0)
    gauss_band.set_matrix(0, 0, 0, 1, -0.5);  /// coefficient of e(0, 1) in the equation for e(0, 0)

    auto res = gauss_band.solve();
    assert(abs(res[0][0] - 4) < 1e-9);  /// it takes 4 moves on average from cell(0, 0)
    assert(abs(res[0][1] - 3) < 1e-9);  /// and 3 moves on average from cell(0, 1)
    assert(abs(res[1][0] - 3) < 1e-9);  /// similarly 3 moves on average from cell(1, 0)
    assert(abs(res[1][1] - 0) < 1e-9);  /// 0 moves for destination

    return 0;
}
