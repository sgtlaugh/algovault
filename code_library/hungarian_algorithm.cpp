/***
 *
 * Hungarian Algorithm (aka Kuhn-Munkres) for weighted bipartite matching or the assignment problem
 * Returns the optimal assignment in a weighted bipartite graph
 * Sum of edges in the optimal matching is either minimum or maximum
 *
 * Let n = number of rows, m = number of columns
 * Complexity: ~O(max(n, m)^3)
 * The actual complexity is O(n^2 * max(n,m))
 *
 * So if n is large and m is small, it would be faster to transpose the matrix before calling the method
 * The difference can be huge:
 *   If n=1000, m=10 then without transpose it would take roughly a second (10^9 operations)
 *   However, with the transpose it becomes O(10^2 * 1000) = 10^5, which will take a few ms
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <class T>
pair<T, vector<pair<int, int>>> hungarian(vector<vector<T>> mat, bool minimize=true){
    if (!mat.size()) return {(T)0, {}};
    int i, j, a, b, d, n, m;

    int row = mat.size(), col = mat[0].size();
    for (i = 0; i < row; i++){
        assert((int)mat[i].size() == col);
    }

    n = row, m = max(n, col);
    vector <int> way(m + 1, 0), match(m + 1, 0), visited(m + 1, 0);
    vector <T> U(m + 1, 0), V(m + 1, 0), P(m + 1, 0), minv(m + 1, 0);

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
                    T v = 0;
                    if (a >= 1 && a <= row && j >= 1 && j <= col) {
                        v = (minimize ? mat[a - 1][j - 1] : -mat[a - 1][j - 1]);
                    }
                    T r = v - U[a] - V[j];

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

    vector <pair<int, int>> matches;

    for (j = 1; j <= m; j++) match[P[j]] = j;
    for (i = 1; i <= m; i++){
        j = match[i];
        if (i >= 1 && i <= row && j >= 1 && j <= col){
            matches.push_back({i - 1, j - 1});
        }
    }
    assert((int)matches.size() == min(row, col));

    T res = minimize ? -V[0] : V[0];
    return {res, matches};
}

int main(){
    int res;
    vector<vector<int>> mat;
    vector<pair<int, int>> matches, exp_matches;

    mat = {
        {1},
    };
    exp_matches = {
        {0, 0},
    };
    tie(res, matches) = hungarian(mat, false);
    assert(res == 1 && matches == exp_matches);

    mat = {
        {1, 2},
        {4, 8},
    };
    exp_matches = {
        {0, 1},
        {1, 0},
    };
    tie(res, matches) = hungarian(mat, true);
    assert(res == 6 && matches == exp_matches);

    mat = {
        {1, 2},
        {4, 8},
    };
    exp_matches = {
        {0, 0},
        {1, 1},
    };
    tie(res, matches) = hungarian(mat, false);
    assert(res == 9 && matches == exp_matches);

    mat = {
        {1, 2},
        {4, 8},
        {16, 32},
        {64, 128},
    };
    exp_matches = {
        {2, 0},
        {3, 1},
    };
    tie(res, matches) = hungarian(mat, false);
    assert(res == 144 && matches == exp_matches);

    mat = {
        {1, 2, 4, 8},
        {16, 32, 48, 64},
    };
    exp_matches = {
        {0, 1}, /// match row 0 with col 1 (2 cost)
        {1, 0}, /// match row 1 with col 0 (16 cost)
    };
    tie(res, matches) = hungarian(mat, true);
    assert(res == 18 && matches == exp_matches);

    mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9},
    };
    exp_matches = {
        {0, 0},
        {1, 1},
        {2, 2},
    };
    tie(res, matches) = hungarian(mat, true);
    assert(res == 15 && matches == exp_matches);

    /// Test with a big matrix to check performance

    int n = 1000;
    mat = vector<vector<int>>(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if ((i + j) % 3 == 0) mat[i][j] = __gcd(i, j);
            else mat[i][j] = i * j;
        }
    }

    clock_t start = clock();
    tie(res, matches) = hungarian(mat, false);
    assert(res == 332833167);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 1.122106 locally
    return 0;
}
