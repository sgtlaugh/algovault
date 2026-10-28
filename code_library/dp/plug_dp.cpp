/***
 *
 * Plug DP
 * Counts or optimizes cycles and paths on a grid: the caller writes one cell's transition, the engine runs the rest
 *
 * Complexity: O(S) expected per cell for S states alive at that cell, plus the transition's own cost
 *
 * PlugDP<T> dp(n, m, bits, start): n x m grid, bits per plug, the all-zero state holds start, equal states add up
 * PlugDP<T, Merge> dp(n, m, bits, start, merge): equal states combine their values with merge(a, b), such as min
 * dp.add(s, v): only inside the transition, offers state s with value v after the current cell
 * dp.run(transition): calls transition(row, col, s, v) for each state s of value v before each cell, row-major, once
 * dp.get(s, pos): plug pos of state s
 * dp.match(s, pos): position of the bracket paired with plug pos, plugs 1 opening and 2 closing, balanced
 * dp.normalize(s): s with its nonzero labels renumbered 1, 2, 3... in order of first appearance from position 0
 * dp.result(missing): value of the all-zero state after the last cell, missing (T() by default) if none
 * dp.set(s, pos, x): state s with plug pos replaced by x, 0 <= x < 2^bits
 * dp.states, dp.values: the states alive before the current cell, or after the last one, and their values
 *
 * Notes:
 *   - Before cell (row, col) plug col is its edge from the left and col + 1 from above, after it down and right
 *   - Plugs are edges: add drops states with an edge leaving the grid, so transitions never test the border
 *   - (m + 1) bits <= 64, bits <= 8: 2 bits fit brackets, labels need 2^bits >= components + 2 to leave a fresh one
 *   - A state holds plugs 0..m only: the row shift moves every bit up by bits, so nothing can be kept above plug m
 *   - main() writes Formula 1 (one cycle, bracket encoding) and HDU 1693 (any cycles, 1 bit) as transitions
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Merge = plus<T>>
struct PlugDP{
    int n, m, bits, row = 0, col = 0, lg = 10;
    uint64_t off_grid = 0;
    Merge merge;
    vector<uint64_t> states, next_states;
    vector<T> values, next_values;
    vector<int> head, link;

    PlugDP(int n, int m, int bits, T start, Merge merge = Merge())
        : n(n), m(m), bits(bits), merge(merge), states{0}, values{start}, head(1 << lg, -1){
        assert(n >= 1 && m >= 1 && bits >= 1 && bits <= 8 && (m + 1) * bits <= 64);
    }

    void add(uint64_t s, T v){
        if (s & off_grid) return;

        size_t h = slot(s);
        for (int i = head[h]; i >= 0; i = link[i]){
            if (next_states[i] == s){
                next_values[i] = merge(next_values[i], v);
                return;
            }
        }

        link.push_back(head[h]);
        head[h] = next_states.size();
        next_states.push_back(s);
        next_values.push_back(v);
    }

    template <typename Transition>
    void run(Transition transition){
        uint64_t ones = (uint64_t(1) << bits) - 1;
        for (row = 0; row < n; row++){
            for (col = 0; col < m; col++){
                /// sized per cell, not per insert: a resize check inside add stops GCC inlining it, 1.6x slower
                while (head.size() < 8 * states.size()) head.assign(size_t(1) << ++lg, -1);
                off_grid = (col == m - 1 ? ones << (m * bits) : 0) | (row == n - 1 ? ones << (col * bits) : 0);
                for (size_t i = 0; i < states.size(); i++) transition(row, col, states[i], values[i]);

                for (uint64_t s : next_states) head[slot(s)] = -1;
                link.clear();
                swap(states, next_states), swap(values, next_values);
                next_states.clear(), next_values.clear();
            }
            for (uint64_t& s : states) s <<= bits;
        }
    }

    int get(uint64_t s, int pos) const{
        return s >> (pos * bits) & ((1 << bits) - 1);
    }

    int match(uint64_t s, int pos) const{
        int dir = get(s, pos) == 1 ? 1 : -1, depth = 0;
        for (int p = pos; p >= 0 && p <= m; p += dir){
            int x = get(s, p);
            depth += (x == 1) - (x == 2);
            if (depth == 0) return p;
        }

        assert(false);
        return -1;
    }

    uint64_t normalize(uint64_t s) const{
        unsigned char id[256] = {};
        int labels = 0;
        for (int p = 0; p <= m; p++){
            int x = get(s, p);
            if (!x) continue;
            if (!id[x]) id[x] = ++labels;
            s = set(s, p, id[x]);
        }

        return s;
    }

    T result(T missing = T()) const{
        return states.empty() ? missing : values[0];
    }

    uint64_t set(uint64_t s, int pos, int x) const{
        int shift = pos * bits;
        return (s & ~(uint64_t((1 << bits) - 1) << shift)) | uint64_t(x) << shift;
    }

    size_t slot(uint64_t s) const{
        return s * 0x9E3779B97F4A7C15ULL >> (64 - lg);
    }
};

int main(){
    /// Formula 1 (Timus 1519): count the single cycles through every free cell, '#' blocked, bracket encoding
    auto formula1 = [](const vector<string>& grid){
        int n = grid.size(), m = grid[0].size(), last = -1;
        for (int i = 0; i < n * m; i++){
            if (grid[i / m][i % m] == '.') last = i;
        }
        if (last == -1) return 0LL;

        PlugDP<long long> dp(n, m, 2, 1);
        dp.run([&](int i, int j, uint64_t s, long long v){
            int l = dp.get(s, j), u = dp.get(s, j + 1);
            uint64_t t = dp.set(dp.set(s, j, 0), j + 1, 0);
            if (grid[i][j] == '#'){
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
            else if (i * m + j == last) dp.add(t, v);  /// closing earlier would leave free cells off the cycle
        });

        return dp.result();
    };
    assert(formula1({"....", "....", "....", "...."}) == 6);
    assert(formula1({"....", "...."}) == 1);
    assert(formula1({"...", ".#.", "..."}) == 1);
    assert(formula1({"...", "..."}) == 1);

    /// HDU 1693: cover every free cell with any number of cycles, one bit per plug since cycles never need pairing
    auto cycle_covers = [](const vector<vector<int>>& open){
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
    };
    assert(cycle_covers({{1, 1, 1, 1}, {1, 1, 1, 1}}) == 2);  /// the ring, or two 2 x 2 squares
    assert(cycle_covers({{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}) == 1);
    assert(cycle_covers({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}) == 0);  /// a cycle on the grid has even length

    return 0;
}
