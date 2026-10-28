#include "../common.h"

#define main library_main
#include "../../code_library/dp/plug_dp.cpp"
#undef main

using Grid = vector<vector<int>>;

long long formula1(const Grid& open){
    int n = open.size(), m = open[0].size(), last = -1;
    for (int i = 0; i < n * m; i++){
        if (open[i / m][i % m]) last = i;
    }
    if (last == -1) return 0;

    PlugDP<long long> dp(n, m, 2, 1);
    dp.run([&](int i, int j, uint64_t s, long long v){
        int l = dp.get(s, j), u = dp.get(s, j + 1);
        uint64_t t = dp.set(dp.set(s, j, 0), j + 1, 0);
        if (!open[i][j]){
            if (!l && !u) dp.add(s, v);
        }
        else if (!l && !u) dp.add(dp.set(dp.set(t, j, 1), j + 1, 2), v);
        else if (!l || !u){
            dp.add(dp.set(t, j, l | u), v);
            dp.add(dp.set(t, j + 1, l | u), v);
        }
        else if (l == 1 && u == 1) dp.add(dp.set(t, dp.match(s, j + 1), 1), v);
        else if (l == 2 && u == 2) dp.add(dp.set(t, dp.match(s, j), 2), v);
        else if (l == 2 && u == 1) dp.add(t, v);
        else if (i * m + j == last) dp.add(t, v);
    });

    /// no edge leaves the bottom row, so the empty profile is the only possible survivor
    assert(dp.states.size() <= 1 && (dp.states.empty() || dp.states[0] == 0));
    return dp.result();
}

long long cycle_covers(const Grid& open){
    int n = open.size(), m = open[0].size();
    PlugDP<long long> dp(n, m, 1, 1);
    dp.run([&](int i, int j, uint64_t s, long long v){
        int l = dp.get(s, j), u = dp.get(s, j + 1);
        uint64_t t = dp.set(dp.set(s, j, 0), j + 1, 0);
        if (!open[i][j]){
            if (!l && !u) dp.add(s, v);
        }
        else if (l && u) dp.add(t, v);
        else if (l || u){
            dp.add(dp.set(t, j, 1), v);
            dp.add(dp.set(t, j + 1, 1), v);
        }
        else dp.add(dp.set(dp.set(t, j, 1), j + 1, 1), v);
    });

    return dp.result();
}

/// Simple paths from cell a to cell b with component labels, each path piece one label; gain(v, i, j) puts (i, j) on it
template <typename T, typename Merge, typename Gain>
T paths(const Grid& open, int a, int b, T start, T missing, Merge merge, Gain gain){
    int n = open.size(), m = open[0].size();
    const int fresh = 15;
    PlugDP<T, Merge> dp(n, m, 4, start, merge);
    dp.run([&](int i, int j, uint64_t s, T v){
        int l = dp.get(s, j), u = dp.get(s, j + 1);
        uint64_t t = dp.set(dp.set(s, j, 0), j + 1, 0);
        bool end = i * m + j == a || i * m + j == b;
        if (!open[i][j]){
            if (!l && !u) dp.add(s, v);
            return;
        }

        if (end){
            if (!l && !u){
                dp.add(dp.normalize(dp.set(t, j, fresh)), gain(v, i, j));
                dp.add(dp.normalize(dp.set(t, j + 1, fresh)), gain(v, i, j));
            }
            else if (!l || !u) dp.add(dp.normalize(t), gain(v, i, j));
        }
        else if (!l && !u){
            dp.add(s, v);
            dp.add(dp.normalize(dp.set(dp.set(t, j, fresh), j + 1, fresh)), gain(v, i, j));
        }
        else if (!l || !u){
            dp.add(dp.set(t, j, l | u), gain(v, i, j));
            dp.add(dp.set(t, j + 1, l | u), gain(v, i, j));
        }
        else if (l != u){
            for (int p = 0; p <= m; p++){
                if (dp.get(t, p) == u) t = dp.set(t, p, l);
            }
            dp.add(dp.normalize(t), gain(v, i, j));
        }
    });

    return dp.result(missing);
}

/// Every edge set where each free cell has degree 2, chosen cell by cell; single keeps only the connected ones
long long brute_cycles(const Grid& open, bool single){
    int n = open.size(), m = open[0].size();
    Grid right(n, vector<int>(m)), down(n, vector<int>(m));
    long long count = 0;

    auto connected = [&](){
        vector<int> parent(n * m);
        iota(parent.begin(), parent.end(), 0);
        function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };
        int cells = 0, parts = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                if (right[i][j]) parent[find(i * m + j)] = find(i * m + j + 1);
                if (down[i][j]) parent[find(i * m + j)] = find((i + 1) * m + j);
            }
        }
        for (int c = 0; c < n * m; c++){
            if (open[c / m][c % m]) cells++, parts += find(c) == c;
        }

        return cells > 0 && parts == 1;
    };

    function<void(int)> go = [&](int c){
        if (c == n * m){
            if (!single || connected()) count++;
            return;
        }

        int i = c / m, j = c % m;
        if (!open[i][j]) return go(c + 1);
        int deg = (j > 0 && right[i][j - 1]) + (i > 0 && down[i - 1][j]);
        bool can_right = j + 1 < m && open[i][j + 1], can_down = i + 1 < n && open[i + 1][j];
        for (int r = 0; r <= can_right; r++){
            for (int d = 0; d <= can_down; d++){
                if (deg + r + d != 2) continue;
                right[i][j] = r, down[i][j] = d;
                go(c + 1);
                right[i][j] = down[i][j] = 0;
            }
        }
    };
    go(0);

    return count;
}

/// Depth-first search over simple paths from a to b: their count and the largest sum of weights on one
pair<long long, long long> brute_paths(const Grid& open, const Grid& weight, int a, int b){
    int n = open.size(), m = open[0].size();
    vector<int> seen(n * m);
    long long count = 0, best = LLONG_MIN;
    function<void(int, long long)> go = [&](int c, long long sum){
        if (c == b){
            count++, best = max(best, sum);
            return;
        }

        int i = c / m, j = c % m, di[] = {0, 0, 1, -1}, dj[] = {1, -1, 0, 0};
        for (int k = 0; k < 4; k++){
            int ni = i + di[k], nj = j + dj[k];
            if (ni < 0 || ni >= n || nj < 0 || nj >= m || !open[ni][nj] || seen[ni * m + nj]) continue;
            seen[ni * m + nj] = 1;
            go(ni * m + nj, sum + weight[ni][nj]);
            seen[ni * m + nj] = 0;
        }
    };
    seen[a] = 1;
    go(a, weight[a / m][a % m]);

    return {count, best};
}

Grid random_grid(int n, int m, int blocked_percent){
    Grid open(n, vector<int>(m));
    for (auto& row : open){
        for (auto& x : row) x = stress::rand_int(1, 100) > blocked_percent;
    }

    return open;
}

int main(){
    /// Hamiltonian cycles of the full 2k x 2k grid, OEIS A003763; 12 x 12 is Timus 1519's largest and needs 60 bits
    vector<long long> full = {1, 6, 1072, 4638576, 467260456608LL, 1076226888605605706LL};
    for (int k = 1; k <= 6; k++) assert(formula1(Grid(2 * k, vector<int>(2 * k, 1))) == full[k - 1]);

    /// Random obstacles on grids up to 6 x 6, single rows and columns included, every free cell on a cycle
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 4 ? 5 : 6), m = stress::rand_int(1, it % 4 ? 5 : 6);
        Grid open = random_grid(n, m, stress::rand_int(0, 40));
        assert(cycle_covers(open) == brute_cycles(open, false));
        assert(formula1(open) == brute_cycles(open, true));
    }

    /// Paths between two free cells: plus counts them, max finds the heaviest, both through component labels
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, 5), m = stress::rand_int(1, it % 3 ? 4 : 6);
        if (n * m > 20) n = 20 / m;
        Grid open = random_grid(n, m, stress::rand_int(0, 35)), weight(n, vector<int>(m));
        for (auto& row : weight){
            for (auto& x : row) x = stress::rand_int(-9, 9);
        }
        vector<int> cells;
        for (int c = 0; c < n * m; c++){
            if (open[c / m][c % m]) cells.push_back(c);
        }
        if (cells.size() < 2) continue;

        int a = cells[stress::rand_int(0, cells.size() - 1)], b = a;
        while (b == a) b = cells[stress::rand_int(0, cells.size() - 1)];
        auto [count, best] = brute_paths(open, weight, a, b);
        auto keep = [](long long v, int, int){ return v; };
        auto add_weight = [&](long long v, int i, int j){ return v + weight[i][j]; };
        auto larger = [](long long x, long long y){ return max(x, y); };
        assert(paths(open, a, b, 1LL, 0LL, plus<long long>(), keep) == count);
        assert(paths(open, a, b, 0LL, LLONG_MIN, larger, add_weight) == best);
    }

    return 0;
}
