/***
 *
 * Coordinate Compression
 * Replaces each value of a vector in place with a small id, and returns the distinct values
 *
 * Complexity: O(n log n), no hashing
 *
 * compress(v): ids follow value order, returns the distinct values sorted, so vals[v[i]] is the original v[i]
 * compress(v, false): ids follow first appearance, returns the distinct values in that order
 *
 * Ids run from 0 to the number of distinct values minus 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
vector<T> compress(vector<T>& v, bool make_sorted=true){
    int n = v.size();
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());

    vector<T> vals;
    vector<int> ids(n);
    for (int i = 0; i < n; i++){
        if (!i || order[i - 1].first < order[i].first) vals.push_back(order[i].first);
        ids[order[i].second] = vals.size() - 1;
    }

    /// Relabel by scanning in input order instead of hashing, std::hash is the identity on ints and hackable
    if (!make_sorted){
        vector<int> relabel(vals.size(), -1);
        vector<T> first_seen;
        for (int& id: ids){
            if (relabel[id] == -1){
                relabel[id] = first_seen.size();
                first_seen.push_back(vals[id]);
            }
            id = relabel[id];
        }
        vals.swap(first_seen);
    }

    for (int i = 0; i < n; i++) v[i] = ids[i];
    return vals;
}

int main(){
    freopen("crazy.in", "r", stdin);
    freopen("crazy.out", "w", stdout);

    int n, c;
    scanf("%d %d", &n, &c);

    vector<int> xs(2 * n + c), ys(2 * n + c);
    for (int i = 0; i < n; i++) scanf("%d %d %d %d", &xs[2 * i], &ys[2 * i], &xs[2 * i + 1], &ys[2 * i + 1]);
    for (int i = 2 * n; i < 2 * n + c; i++) scanf("%d %d", &xs[i], &ys[i]);

    int w = 2 * compress(xs).size() + 1, h = 2 * compress(ys).size() + 1;
    enum {FREE, FENCE, SEEN, COW};
    vector<char> grid(w * h, FREE);

    for (int i = 0; i < n; i++){
        int x1 = 2 * xs[2 * i] + 1, y1 = 2 * ys[2 * i] + 1, x2 = 2 * xs[2 * i + 1] + 1, y2 = 2 * ys[2 * i + 1] + 1;
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);
        for (int x = x1; x <= x2; x++){
            for (int y = y1; y <= y2; y++) grid[x * h + y] = FENCE;
        }
    }
    for (int i = 2 * n; i < 2 * n + c; i++) grid[(2 * xs[i] + 1) * h + 2 * ys[i] + 1] = COW;

    int res = 0;
    vector<int> stk;
    for (int i = 2 * n; i < 2 * n + c; i++){
        int s = (2 * xs[i] + 1) * h + 2 * ys[i] + 1;
        if (grid[s] != COW) continue;

        int cows = 1;
        grid[s] = SEEN;
        stk.assign(1, s);
        while (!stk.empty()){
            int u = stk.back(), x = u / h, y = u % h;
            stk.pop_back();

            int dx[] = {1, -1, 0, 0}, dy[] = {0, 0, 1, -1};
            for (int d = 0; d < 4; d++){
                int nx = x + dx[d], ny = y + dy[d], v = nx * h + ny;
                if (nx < 0 || nx >= w || ny < 0 || ny >= h || grid[v] == FENCE || grid[v] == SEEN) continue;
                cows += grid[v] == COW;
                grid[v] = SEEN;
                stk.push_back(v);
            }
        }
        res = max(res, cows);
    }

    printf("%d\n", res);
    return 0;
}
