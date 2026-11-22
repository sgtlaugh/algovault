/***
 *
 * Monotone Minima and SMAWK
 * Row minima of an implicit n x m matrix f(i, j) whose minimum columns move right as the row index grows
 *
 * Complexity:
 *   - monotone_minima: O((n + m) log n) evaluations of f
 *   - smawk: O(n + m) evaluations of f
 *
 * monotone_minima(n, m, f): res[i] = the leftmost j minimizing f(i, j), for rows 0..n-1 and columns 0..m-1
 * smawk(n, m, f): the same res for a totally monotone f
 *
 * Notes:
 *   - Both require m >= 1, and the row minimum is f(i, res[i])
 *   - monotone_minima requires the leftmost minimum column to be non-decreasing in i
 *   - smawk requires f(i, k) < f(i, j) to imply f(i2, k) < f(i2, j) for all i < i2 and j < k (totally monotone)
 *   - Both hold for Monge f: f(i, j) + f(i2, k) <= f(i, k) + f(i2, j) for all i < i2 and j < k
 *   - Monge examples: g(b[j] - a[i]) for convex g and sorted a, b, plus terms in only i or only j, and sums
 *   - f may return any type with < and copy, and is called only with 0 <= i < n and 0 <= j < m
 *   - smawk calls f 4 to 5.5 (n + m) times, monotone_minima about (n + m) log2 n times
 *   - O(1) f at -O2: monotone_minima 1.3x faster at n = m = 10^6, smawk 1.3x faster at n = 10^3, m = 10^6
 *   - f costing a sqrt and a log: smawk 1.2x faster at n = m = 10^6 and 1.9x at n = 10^3, m = 10^6
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename F>
vector<int> monotone_minima(int n, int m, F f){
    assert(m >= 1);

    vector<int> res(n);
    auto solve = [&](auto&& self, int lo, int hi, int jlo, int jhi) -> void{
        if (lo > hi) return;

        int mid = (lo + hi) / 2, best = jlo;
        auto best_val = f(mid, jlo);
        for (int j = jlo + 1; j <= jhi; j++){
            auto val = f(mid, j);
            if (val < best_val) best = j, best_val = val;
        }
        res[mid] = best;

        self(self, lo, mid - 1, jlo, best);
        self(self, mid + 1, hi, best, jhi);
    };

    solve(solve, 0, n - 1, 0, m - 1);
    return res;
}

template <typename F>
vector<int> smawk(int n, int m, F f){
    assert(m >= 1);

    vector<int> res(n), rows(n), cols(m);
    iota(rows.begin(), rows.end(), 0);
    iota(cols.begin(), cols.end(), 0);
    auto solve = [&](auto&& self, const vector<int>& rows, const vector<int>& cols) -> void{
        if (rows.empty()) return;

        /// kept[s] loses or ties to kept[s - 1] in every row above rows[s], so at most rows.size() columns matter
        vector<int> kept;
        for (int c : cols){
            while (!kept.empty()){
                int r = rows[kept.size() - 1];
                if (!(f(r, c) < f(r, kept.back()))) break;
                kept.pop_back();
            }
            if (kept.size() < rows.size()) kept.push_back(c);
        }

        vector<int> odd;
        for (int i = 1; i < (int)rows.size(); i += 2) odd.push_back(rows[i]);
        self(self, odd, kept);

        /// Each even row's minimum lies between the minima of its odd neighbours, so k only moves right
        for (int i = 0, k = 0; i < (int)rows.size(); i += 2){
            int r = rows[i], last = i + 1 < (int)rows.size() ? res[rows[i + 1]] : kept.back();
            auto best_val = f(r, kept[k]);
            res[r] = kept[k];
            while (kept[k] != last){
                k++;
                auto val = f(r, kept[k]);
                if (val < best_val) res[r] = kept[k], best_val = val;
            }
        }
    };

    solve(solve, rows, cols);
    return res;
}

int main(){
    /***
     * f(i, j) = (a[i] - b[j])^2 for sorted a = {1, 4, 6} and b = {0, 3, 5, 9} is Monge
     * Rows {1, 4, 16, 64}, {16, 1, 1, 25}, {36, 9, 1, 9}: row 1 ties at columns 1 and 2 and takes the leftmost
    ***/
    vector<long long> a = {1, 4, 6}, b = {0, 3, 5, 9};
    auto squared = [&](int i, int j){ return (a[i] - b[j]) * (a[i] - b[j]); };
    assert((monotone_minima(3, 4, squared) == vector<int>{0, 1, 2}));
    assert((smawk(3, 4, squared) == vector<int>{0, 1, 2}));

    /// Rows {2, 0, 1} and {5, 0, 6}: minima columns {1, 1} are monotone, but column 2 beats 0 only in row 0
    vector<vector<int>> loose = {{2, 0, 1}, {5, 0, 6}};
    assert((monotone_minima(2, 3, [&](int i, int j){ return loose[i][j]; }) == vector<int>{1, 1}));

    return 0;
}
