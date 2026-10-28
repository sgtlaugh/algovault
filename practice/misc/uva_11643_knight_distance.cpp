/***
 *
 * Knight Distance
 * Fewest knight moves between two squares of an infinite chessboard (every quadrant usable)
 *
 * Complexity: O(1)
 *
 * knight_distance(x, y): moves from (0, 0) to (x, y), any long long coordinates with |x|, |y| <= 1e18
 * For two squares use knight_distance(x2 - x1, y2 - y1)
 *
 * By symmetry only x >= y >= 0 matters. Far from the diagonal (2y < x) the knight walks along x with
 * (2, +-1) steps, near the diagonal it mixes (2, 1) and (1, 2) steps. (1, 0) and (2, 2) are the two exceptions
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long knight_distance(long long x, long long y){
    x = llabs(x), y = llabs(y);
    if (x < y) swap(x, y);
    if (x == 1 && y == 0) return 3;
    if (x == 2 && y == 2) return 4;

    if (y == 0 || 2 * y < x){
        long long c = y & 1, a = x - 2 * c, b = a & 3;
        return (a - b) / 2 + b + c;
    }

    long long d = x - (x - y) / 2, c = (x - y) & 1, z = d % 3 != 0;
    return d / 3 * 2 + c + z * 2 * (1 - c);
}

int main(){
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, k;
        scanf("%d %d", &n, &k);
        vector<int> r(k), c(k);
        for (int i = 0; i < k; i++) scanf("%d %d", &r[i], &c[i]);

        auto corner = [&](int i){ return (r[i] == 1 || r[i] == n) && (c[i] == 1 || c[i] == n); };
        vector<vector<int>> dist(k, vector<int>(k));
        for (int i = 0; i < k; i++){
            for (int j = 0; j < k; j++){
                dist[i][j] = knight_distance(r[i] - r[j], c[i] - c[j]);
                if (abs(r[i] - r[j]) == 1 && abs(c[i] - c[j]) == 1 && (corner(i) || corner(j))) dist[i][j] = 4;
                if (n == 4 && corner(i) && corner(j) && (r[i] == r[j]) != (c[i] == c[j])) dist[i][j] = 5;
            }
        }

        int full = 1 << k;
        vector<int> dp(full * k, INT_MAX / 2);
        dp[k] = 0;
        for (int mask = 3; mask < full; mask += 2){
            for (int rest = mask & ~1; rest; rest &= rest - 1){
                int j = __builtin_ctz(rest), prev = (mask ^ 1 << j) * k, best = INT_MAX / 2;
                for (int from = mask ^ 1 << j; from; from &= from - 1){
                    int i = __builtin_ctz(from);
                    best = min(best, dp[prev + i] + dist[i][j]);
                }
                dp[mask * k + j] = best;
            }
        }

        int best = k == 1 ? 0 : INT_MAX;
        for (int i = 1; i < k; i++) best = min(best, dp[(full - 1) * k + i] + dist[i][0]);
        printf("Case %d: %d\n", tc, best);
    }

    return 0;
}
