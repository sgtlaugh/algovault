/***
 *
 * Maximum Square and Diamond
 * Largest all-filled square and diamond ending at every cell of a 0/1 grid
 *
 * Complexity: O(rows * cols)
 *
 * square_sizes(g)[i][j]: side of the largest filled square whose bottom-right cell is (i, j)
 * diamond_sizes(g)[i][j]: size of the largest filled diamond whose bottom tip is (i, j)
 *     a diamond of size k is every cell within Manhattan distance k - 1 of (i - k + 1, j)
 * A value k counts k shapes ending there (sizes 1..k), so the sum is the number of filled squares or diamonds
 * and the maximum is the largest one; g[i][j] != 0 means filled
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<vector<int>> square_sizes(const vector<vector<int>>& g){
    int n = g.size(), m = n ? g[0].size() : 0;
    vector<vector<int>> dp(n, vector<int>(m, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if (!g[i][j]) continue;
            int up = i ? dp[i - 1][j] : 0, left = j ? dp[i][j - 1] : 0, x = min(up, left);
            dp[i][j] = x + (g[i - x][j - x] != 0);  /// the corner cell decides between min(up, left) and one more
        }
    }
    return dp;
}

vector<vector<int>> diamond_sizes(const vector<vector<int>>& g){
    int n = g.size(), m = n ? g[0].size() : 0;
    vector<vector<int>> dp(n, vector<int>(m, 0));
    auto filled = [&](int i, int j){ return i >= 0 && j >= 0 && j < m && g[i][j]; };
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if (!g[i][j]) continue;
            int a = (i && j) ? dp[i - 1][j - 1] : 0, b = (i && j + 1 < m) ? dp[i - 1][j + 1] : 0, x = min(a, b);
            /// Two diamonds of size x on the upper diagonals cover all but the two cells on the axis above them
            if (!x || !filled(i - 1, j)) dp[i][j] = 1;
            else if (filled(i - 2 * x, j) && filled(i - 2 * x + 1, j)) dp[i][j] = x + 1;
            else dp[i][j] = x;
        }
    }
    return dp;
}

int main(){
    vector<vector<int>> g = {
        {1, 1, 1, 0},
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {0, 1, 1, 1},
    };
    auto sq = square_sizes(g);
    assert(sq[2][2] == 3 && sq[3][3] == 3 && sq[3][0] == 0 && sq[0][0] == 1 && sq[1][1] == 2 && sq[1][3] == 1);

    vector<vector<int>> d = {
        {0, 0, 1, 0, 0},
        {0, 1, 1, 1, 0},
        {1, 1, 1, 1, 1},
        {0, 1, 1, 1, 0},
        {0, 0, 1, 0, 0},
    };
    auto dm = diamond_sizes(d);
    assert(dm[4][2] == 3 && dm[3][2] == 2 && dm[2][2] == 2 && dm[0][2] == 1 && dm[2][0] == 1);

    assert(square_sizes({}).empty() && diamond_sizes({}).empty());
    return 0;
}
