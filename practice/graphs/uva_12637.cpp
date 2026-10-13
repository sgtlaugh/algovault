/***
 *
 * Hungarian Algorithm (aka Kuhn-Munkres) for weighted bipartite matching or the assignment problem
 * Matches rows to distinct columns of a cost matrix so the sum of the chosen cells is minimum or maximum
 *
 * auto [total, pairs] = hungarian(mat, minimize = true): mat is an n x m vector<vector<T>>, rectangular allowed
 *     pairs holds min(n, m) (row, col) pairs, 0-based, sorted by row, each row and column appears at most once
 *     total is the sum of mat[row][col] over pairs, minimize = false maximizes it instead
 *     T must be signed (int, long long or double), 2 * max(n, m) * max|mat[i][j]| must fit in T
 *     (potentials reach about max(n, m) * (max - min) of the entries, 4C on a 3 x 3 matrix of +-C)
 *     An empty mat returns {0, {}}
 *
 * Let n = number of rows, m = number of columns
 * Complexity: O(min(n, m)^2 * max(n, m))
 * Tall matrices (n > m) are transposed internally, so the smaller side is always the one iterated
 * With ties, a tall matrix may get a different optimal matching than its transpose would, the cost is the same
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
pair<T, vector<pair<int, int>>> hungarian(vector<vector<T>> mat, bool minimize=true){
    if (!mat.size()) return {(T)0, {}};
    int i, j, a, b, d, n, m;

    n = mat.size(), m = mat[0].size();
    for (i = 0; i < n; i++){
        assert((int)mat[i].size() == m);
    }

    if (n > m){
        vector<vector<T>> t(m, vector<T>(n));
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++) t[j][i] = mat[i][j];
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

    vector<int> way(m + 1, 0), match(n + 1, 0), visited(m + 1, 0), P(m + 1, 0);
    vector<T> U(n + 1, 0), V(m + 1, 0), minv(m + 1, 0);

    for (i = 1; i <= n; i++){
        b = 0, P[0] = i;
        for (j = 0; j <= m; j++){
            minv[j] = numeric_limits<T>::max(), visited[j] = 0;
        }
        do{
            T w = numeric_limits<T>::max();
            a = P[b], d = 0, visited[b] = 1;

            for (j = 1; j <= m; j++){
                if (!visited[j]){
                    T r = mat[a - 1][j - 1] - U[a] - V[j];

                    if (r < minv[j]) minv[j] = r, way[j] = b;
                    if (minv[j] < w) w = minv[j], d = j;
                }
            }

            for (j = 0, b = d; j <= m; j++){
                if (visited[j]) U[P[j]] += w, V[j] -= w;
                else minv[j] -= w;
            }
        } while (P[b]);

        do{
            d = way[b];
            P[b] = P[d], b = d;
        } while (b);
    }

    vector<pair<int, int>> matches;

    for (j = 1; j <= m; j++) match[P[j]] = j;
    for (i = 1; i <= n; i++) matches.push_back({i - 1, match[i] - 1});

    return {minimize ? -V[0] : V[0], matches};
}

int main(){
    int r, n;
    while (scanf("%d %d", &r, &n) == 2){
        vector<int> x(r + n), y(r + n);
        for (int i = 0; i < r + n; i++) scanf("%d %d", &x[i], &y[i]);

        /// Each delivery guy makes one round from his own restaurant, so orders go to distinct restaurants
        vector<vector<int>> dist(r, vector<int>(n));
        for (int i = 0; i < r; i++){
            for (int j = 0; j < n; j++) dist[i][j] = abs(x[i] - x[r + j]) + abs(y[i] - y[r + j]);
        }
        printf("%d\n", hungarian(dist).first);
    }
    return 0;
}
