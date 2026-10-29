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
        int l = dp.left(s), u = dp.up(s);
        if (!open[i][j]){
            if (!l && !u) dp.add(s, v);
        }
        else if (!l && !u) dp.add(dp.place(s, 1, 2), v);
        else if (!l || !u){
            dp.add(dp.place(s, l | u, 0), v);
            dp.add(dp.place(s, 0, l | u), v);
        }
        else if (l == 1 && u == 2){
            if (i * m + j == last) dp.add(dp.place(s, 0, 0), v);
        }
        else dp.add(dp.join(s), v);
    });

    /// no edge leaves the bottom row, so the empty profile is the only possible survivor
    assert(dp.states.size() <= 1 && (dp.states.empty() || dp.states[0] == 0));
    return dp.result();
}

/// Raw get and set on purpose, so transitions written without left, up and place stay covered
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
T paths(const Grid& open, int a, int b, T start, T missing, Merge merge, Gain gain, int bits = 4){
    /// The top bit as a fresh label: above any label these grids reach, and lost by code that keeps only 4 bits
    int n = open.size(), m = open[0].size(), fresh = 1 << (bits - 1);
    PlugDP<T, Merge> dp(n, m, bits, start, merge);
    dp.run([&](int i, int j, uint64_t s, T v){
        int l = dp.left(s), u = dp.up(s);
        bool end = i * m + j == a || i * m + j == b;
        if (!open[i][j]){
            if (!l && !u) dp.add(s, v);
            return;
        }

        if (end){
            if (!l && !u){
                dp.add(dp.normalize(dp.place(s, fresh, 0)), gain(v, i, j));
                dp.add(dp.normalize(dp.place(s, 0, fresh)), gain(v, i, j));
            }
            else if (!l || !u) dp.add(dp.normalize(dp.place(s, 0, 0)), gain(v, i, j));
        }
        else if (!l && !u){
            dp.add(s, v);
            dp.add(dp.normalize(dp.place(s, fresh, fresh)), gain(v, i, j));
        }
        else if (!l || !u){
            dp.add(dp.place(s, l | u, 0), gain(v, i, j));
            dp.add(dp.place(s, 0, l | u), gain(v, i, j));
        }
        else if (l != u){
            uint64_t t = dp.place(s, 0, 0);
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

/// Plus counts the paths between two random free cells, max finds the heaviest, both through component labels
void check_paths(const Grid& open, int bits){
    int n = open.size(), m = open[0].size();
    Grid weight(n, vector<int>(m));
    for (auto& row : weight){
        for (auto& x : row) x = stress::rand_int(-9, 9);
    }
    vector<int> cells;
    for (int c = 0; c < n * m; c++){
        if (open[c / m][c % m]) cells.push_back(c);
    }
    if (cells.size() < 2) return;

    int a = cells[stress::rand_int(0, cells.size() - 1)], b = a;
    while (b == a) b = cells[stress::rand_int(0, cells.size() - 1)];
    auto [count, best] = brute_paths(open, weight, a, b);
    auto keep = [](long long v, int, int){ return v; };
    auto add_weight = [&](long long v, int i, int j){ return v + weight[i][j]; };
    auto larger = [](long long x, long long y){ return max(x, y); };
    assert(paths(open, a, b, 1LL, 0LL, plus<long long>(), keep, bits) == count);
    assert(paths(open, a, b, 0LL, LLONG_MIN, larger, add_weight, bits) == best);
}

/// Pairs of a bracket sequence found with a stack, asserting the sequence is balanced
set<pair<int, int>> bracket_pairs(const vector<int>& seq){
    vector<int> opened;
    set<pair<int, int>> pairs;
    for (int p = 0; p < (int)seq.size(); p++){
        assert(seq[p] >= 0 && seq[p] <= 2);
        if (seq[p] == 1) opened.push_back(p);
        if (seq[p] == 2){
            assert(!opened.empty());
            pairs.insert({opened.back(), p});
            opened.pop_back();
        }
    }
    assert(opened.empty());

    return pairs;
}

/// join at every column it allows: the two entering pairs become one between their far ends, other pairs stay
void check_join(const vector<int>& seq, int bits){
    int m = seq.size() - 1, top = (1 << bits) - 1;
    PlugDP<long long> dp(1, m, bits, 1);
    uint64_t s = 0;
    for (int p = 0; p <= m; p++) s |= uint64_t(seq[p]) << (p * bits);
    set<pair<int, int>> pairs = bracket_pairs(seq);

    for (int col = 0; col < m; col++){
        int l = seq[col], u = seq[col + 1];
        if (!l || !u || (l == 1 && u == 2)) continue;

        set<pair<int, int>> want = pairs;
        vector<int> ends;
        for (auto [a, b] : pairs){
            if (a == col || a == col + 1) ends.push_back(b), want.erase({a, b});
            if (b == col || b == col + 1) ends.push_back(a), want.erase({a, b});
        }
        assert(ends.size() == 2);
        want.insert({min(ends[0], ends[1]), max(ends[0], ends[1])});

        dp.col = col;
        uint64_t joined = dp.join(s);
        vector<int> got(m + 1);
        for (int p = 0; p <= m; p++) got[p] = joined >> (p * bits) & top;
        if ((m + 1) * bits < 64) assert(joined >> ((m + 1) * bits) == 0);
        assert(got[col] == 0 && got[col + 1] == 0);
        assert(bracket_pairs(got) == want);
    }
}

/// left, up and place against get and set at every column, plug values up to 2^bits - 1
void check_cell_calls(int m, int bits){
    PlugDP<long long> dp(1, m, bits, 1);
    int top = (1 << bits) - 1;
    uint64_t s = 0;
    for (int p = 0; p <= m; p++) s = dp.set(s, p, stress::rand_int(0, top));

    for (int col = 0; col < m; col++){
        dp.col = col;
        int down = stress::rand_int(0, top), right = stress::rand_int(0, top);
        assert(dp.left(s) == dp.get(s, col));
        assert(dp.up(s) == dp.get(s, col + 1));
        assert(dp.place(s, down, right) == dp.set(dp.set(s, col, down), col + 1, right));
    }
}

/// join must assert, not return a state, on left 1 with up 2, on a missing plug and on a label that is no bracket
void check_rejects(int l, int u){
    stress::expect_assert("dp/plug_dp.cpp", [&]{
        PlugDP<long long> dp(1, 3, 2, 1);
        dp.col = 1;
        dp.join(dp.set(dp.set(dp.set(dp.set(0, 0, 1), 1, l), 2, u), 3, 2));
    });
}

/// A uniformly chosen next symbol among those that can still end balanced
vector<int> random_brackets(int len){
    vector<int> seq(len);
    int depth = 0;
    for (int p = 0; p < len; p++){
        int remaining = len - p;
        vector<int> options;
        if (depth < remaining) options.push_back(0);
        if (depth + 1 < remaining) options.push_back(1);
        if (depth > 0) options.push_back(2);
        seq[p] = options[stress::rand_int(0, options.size() - 1)];
        depth += (seq[p] == 1) - (seq[p] == 2);
    }

    return seq;
}

Grid random_grid(int n, int m, int blocked_percent){
    Grid open(n, vector<int>(m));
    for (auto& row : open){
        for (auto& x : row) x = stress::rand_int(1, 100) > blocked_percent;
    }

    return open;
}

/// Free cells only in the last w columns, where plug m sits in the top bits once (m + 1) bits = 64
Grid right_window(int n, int m, int w, int blocked_percent){
    Grid open = random_grid(n, m, blocked_percent);
    for (auto& row : open) fill(row.begin(), row.end() - w, 0);

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

    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, 5), m = stress::rand_int(1, it % 3 ? 4 : 6);
        if (n * m > 20) n = 20 / m;
        check_paths(random_grid(n, m, stress::rand_int(0, 35)), 4);
    }

    /// Every bit width at (m + 1) bits = 64: a shift or mask past bit 63, or a get or set cut to 4 bits, breaks these
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 5), w = stress::rand_int(1, 4), blocked = stress::rand_int(0, 35);
        Grid open = right_window(n, 63, w, blocked);
        assert(cycle_covers(open) == brute_cycles(open, false));
        open = right_window(n, 31, w, blocked);
        assert(formula1(open) == brute_cycles(open, true));
        check_paths(right_window(n, 15, w, blocked), 4);
        check_paths(right_window(n, 7, w, blocked), 8);
    }

    for (auto [l, u] : vector<pair<int, int>>{{1, 2}, {0, 0}, {0, 1}, {2, 0}, {3, 1}, {2, 3}}) check_rejects(l, u);

    /// Every balanced sequence of 2 to 8 plugs, at 2 bits and at 3, where a label read past bit 1 shows up
    for (int len = 2; len <= 8; len++){
        int total = 1;
        for (int p = 0; p < len; p++) total *= 3;
        for (int code = 0; code < total; code++){
            vector<int> seq(len);
            int rest = code, depth = 0;
            bool balanced = true;
            for (int& x : seq){
                x = rest % 3;
                rest /= 3;
                depth += (x == 1) - (x == 2);
                balanced &= depth >= 0;
            }
            if (!balanced || depth) continue;
            check_join(seq, 2);
            check_join(seq, 3);
        }
    }

    /// Up to 64 bits of state: plug m in the top bits once (m + 1) bits = 64, and long scans for match
    for (long long it = 0; it < stress::scaled(3000); it++){
        int bits = stress::rand_int(2, 8), m = it % 4 ? stress::rand_int(1, 64 / bits - 1) : 64 / bits - 1;
        check_join(random_brackets(m + 1), bits);
        bits = stress::rand_int(1, 8), m = it % 4 ? stress::rand_int(1, 64 / bits - 1) : 64 / bits - 1;
        check_cell_calls(m, bits);
    }

    return 0;
}
