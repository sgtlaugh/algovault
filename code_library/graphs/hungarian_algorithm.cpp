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

    /// Column potentials start at 0 and only fall, so an unmatched column stays optimal when n < m
    vector<int> col_of_row(n, -1), row_of_col(m, -1), order(m), prev_row(m);
    vector<T> col_potential(m, 0), dist(m);

    for (int start = 0; start < n; start++){
        int best = 0;
        for (int j = 0; j < m; j++){
            dist[j] = mat[start][j] - col_potential[j], prev_row[j] = start, order[j] = j;
            if (dist[j] < dist[best]) best = j;
        }

        /// order[0, done) holds the settled columns, best is the position of the closest unsettled one
        int done = 0, col;
        while (true){
            swap(order[done], order[best]);
            col = order[done++];
            if (row_of_col[col] == -1) break;

            int row = row_of_col[col];
            T base = dist[col] - (mat[row][col] - col_potential[col]), low = numeric_limits<T>::max();
            for (int k = done; k < m; k++){
                int j = order[k];
                T r = mat[row][j] - col_potential[j] + base;
                if (r < dist[j]) dist[j] = r, prev_row[j] = row;

                /// A free column wins a tie, many-tie inputs would otherwise settle long runs of equal columns first
                if (dist[j] < low) low = dist[j], best = k;
                else if (dist[j] == low && row_of_col[j] == -1) best = k;
            }
        }

        for (int k = 0; k < done; k++) col_potential[order[k]] += dist[order[k]] - dist[col];
        for (int row = -1; row != start;){
            row = prev_row[col];
            row_of_col[col] = row;
            swap(col_of_row[row], col);
        }
    }

    T total = 0;
    vector<pair<int, int>> matches;
    for (int i = 0; i < n; i++){
        matches.push_back({i, col_of_row[i]});
        total += mat[i][col_of_row[i]];
    }

    return {minimize ? total : -total, matches};
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
