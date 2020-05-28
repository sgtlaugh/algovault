/***
 *
 * Hungarian Algorithm (aka Kuhn-Munkres) for weighted bipartite matching or the assignment problem
 * Returns the optimal assignment in a weighted bipartite graph where the sum of edges in the matching is either minimum or maximum
 *
 * Complexity: O(n^3), where n = max(row, col)
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
pair<T, vector<pair<int, int>>> hungarian_method(vector<vector<T>> mat, bool minimize=true){
    if (!mat.size()) return {(T)0, {}};
    int i, j, a, b, d, n, m, row = mat.size(), col = mat[0].size();
    for (i = 0; i < row; i++) assert((int)mat[i].size() == col);

    n = row, m = max(n, col);
    vector <int> way(m + 1, 0), match(m + 1, 0), visited(m + 1, 0);
    vector <T> U(m + 1, 0), V(m + 1, 0), P(m + 1, 0), minv(m + 1, 0);

    for (i = 1; i <= n; i++){
        b = 0, P[0] = i;
        for (j = 0; j <= m; j++) minv[j] = numeric_limits<T>::max(), visited[j] = 0;
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

    T res = minimize ? -V[0] : V[0];

    vector <pair<int, int>> matched_pairs;
    for (j = 1; j <= m; j++) match[P[j]] = j;
    for (i = 1; i <= m; i++){
        j = match[i];
        if (i >= 1 && i <= row && j >= 1 && j <= col){
            matched_pairs.push_back({i - 1, j - 1});
        }
    }
    assert((int)matched_pairs.size() == min(row, col));

    return {res, matched_pairs};
}

int main(){
    int res;
    vector<vector<int>> mat;
    vector<pair<int, int>> matched_pairs, expected_matched_pairs;

    mat = {{1}};
    expected_matched_pairs = {{0, 0}};
    tie(res, matched_pairs) = hungarian_method(mat, false);
    assert(res == 1 && matched_pairs == expected_matched_pairs);

    mat = {{1, 2}, {4, 8}};
    expected_matched_pairs = {{0, 1}, {1, 0}};
    tie(res, matched_pairs) = hungarian_method(mat, true);
    assert(res == 6 && matched_pairs == expected_matched_pairs);

    mat = {{1, 2}, {4, 8}};
    expected_matched_pairs = {{0, 0}, {1, 1}};
    tie(res, matched_pairs) = hungarian_method(mat, false);
    assert(res == 9 && matched_pairs == expected_matched_pairs);

    mat = {{1, 2}, {4, 8}, {16, 32}, {64, 128}};
    expected_matched_pairs = {{2, 0}, {3, 1}};
    tie(res, matched_pairs) = hungarian_method(mat, false);
    assert(res == 144 && matched_pairs == expected_matched_pairs);

    mat = {{1, 2, 4, 8}, {16, 32, 48, 64}};
    expected_matched_pairs = {{0, 1}, {1, 0}};
    tie(res, matched_pairs) = hungarian_method(mat, true);
    assert(res == 18 && matched_pairs == expected_matched_pairs);

    mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    expected_matched_pairs = {{0, 0}, {1, 1}, {2, 2}};
    tie(res, matched_pairs) = hungarian_method(mat, true);
    assert(res == 15 && matched_pairs == expected_matched_pairs);

    int n = 1000;
    mat = vector<vector<int>>(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if ((i + j) % 3 == 0) mat[i][j] = __gcd(i, j);
            else mat[i][j] = i * j;
        }
    }

    clock_t start = clock();
    tie(res, matched_pairs) = hungarian_method(mat, false);
    assert(res == 332833167);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 1.122106 locally
    return 0;
}
