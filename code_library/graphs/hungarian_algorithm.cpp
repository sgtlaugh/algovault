/***
 *
 * Hungarian Algorithm (Kuhn-Munkres)
 * Matches rows to distinct columns of a cost matrix so the sum of the chosen cells is minimum or maximum
 *
 * Complexity: O(min(n, m)^2 max(n, m)) for n rows and m columns
 *
 * hungarian(mat, minimize = true): {total, pairs}, an optimal assignment of the n x m vector<vector<T>> mat
 *
 * Notes:
 *   - pairs has min(n, m) (row, col) pairs sorted by row, no row or column twice, total sums their cells
 *   - minimize = false maximizes the total instead
 *   - T is signed (int, long long or double) and 2 max(n, m) max |mat[i][j]| must fit in T
 *   - An empty mat returns {0, {}}
 *   - With ties, a tall matrix may get a different optimal matching than its transpose, the total is the same
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
pair<T, vector<pair<int, int>>> hungarian(vector<vector<T>> mat, bool minimize=true){
    if (!mat.size()) return {(T)0, {}};

    int n = mat.size(), m = mat[0].size();
    for (int i = 0; i < n; i++){
        assert((int)mat[i].size() == m);
    }

    if (n > m){
        vector<vector<T>> t(m, vector<T>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) t[j][i] = mat[i][j];
        }
        auto res = hungarian(t, minimize);
        for (auto& [x, y]: res.second) swap(x, y);
        sort(res.second.begin(), res.second.end());
        return res;
    }

    if (!minimize){
        for (auto& r: mat){
            for (auto& x: r) x = -x;
        }
    }

    vector<int> way(m + 1, 0), match(n + 1, 0), visited(m + 1, 0), row_of_col(m + 1, 0);
    vector<T> row_potential(n + 1, 0), col_potential(m + 1, 0), minv(m + 1, 0);

    for (int i = 1; i <= n; i++){
        int col = 0;
        row_of_col[0] = i;
        for (int j = 0; j <= m; j++){
            minv[j] = numeric_limits<T>::max(), visited[j] = 0;
        }
        do{
            T w = numeric_limits<T>::max();
            int row = row_of_col[col], next = 0;
            visited[col] = 1;

            for (int j = 1; j <= m; j++){
                if (!visited[j]){
                    T r = mat[row - 1][j - 1] - row_potential[row] - col_potential[j];

                    if (r < minv[j]) minv[j] = r, way[j] = col;
                    if (minv[j] < w) w = minv[j], next = j;
                }
            }

            for (int j = 0; j <= m; j++){
                if (visited[j]) row_potential[row_of_col[j]] += w, col_potential[j] -= w;
                else minv[j] -= w;
            }
            col = next;
        } while (row_of_col[col]);

        do{
            int prev = way[col];
            row_of_col[col] = row_of_col[prev], col = prev;
        } while (col);
    }

    vector<pair<int, int>> matches;

    for (int j = 1; j <= m; j++) match[row_of_col[j]] = j;
    for (int i = 1; i <= n; i++) matches.push_back({i - 1, match[i] - 1});

    return {minimize ? -col_potential[0] : col_potential[0], matches};
}

int main(){
    /// Cheapest is 2 + 4 = 6 on the anti-diagonal, dearest 1 + 8 = 9 on the diagonal
    vector<vector<int>> mat = {
        {1, 2},
        {4, 8},
    };
    auto [cheapest, cheap_pairs] = hungarian(mat);
    assert(cheapest == 6);
    assert((cheap_pairs == vector<pair<int, int>>{{0, 1}, {1, 0}}));

    auto [dearest, dear_pairs] = hungarian(mat, false);
    assert(dearest == 9);
    assert((dear_pairs == vector<pair<int, int>>{{0, 0}, {1, 1}}));

    /// A tall matrix matches only 2 of its 4 rows: 16 + 128 beats every other pair of distinct columns
    vector<vector<long long>> tall = {
        {1, 2},
        {4, 8},
        {16, 32},
        {64, 128},
    };
    auto [best, pairs] = hungarian(tall, false);
    assert(best == 144);
    assert((pairs == vector<pair<int, int>>{{2, 0}, {3, 1}}));

    return 0;
}
