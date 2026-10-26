#include "../common.h"

#define main library_main
#include "../../code_library/dp/blocks_dp.cpp"
#undef main

/// Every click order: clicking removes one maximal run of equal colors
int brute(const vector<int>& v, map<vector<int>, int>& memo){
    if (v.empty()) return 0;
    auto found = memo.find(v);
    if (found != memo.end()) return found->second;

    int best = 0;
    for (size_t i = 0, j; i < v.size(); i = j){
        for (j = i; j < v.size() && v[j] == v[i]; j++) {}
        vector<int> rest(v.begin(), v.begin() + i);
        rest.insert(rest.end(), v.begin() + j, v.end());
        best = max(best, (int)((j - i) * (j - i)) + brute(rest, memo));
    }

    return memo[v] = best;
}

/// The previous top-down memo, carry index never capped, so it checks the per-run cap at sizes brute cannot reach
int memo_reference(const vector<int>& v){
    int n = v.size();
    vector<int> color, len;
    for (int i = 0; i < n; i++){
        if (!len.empty() && color.back() == v[i]) len.back()++;
        else color.push_back(v[i]), len.push_back(1);
    }

    int m = len.size();
    vector<int> dp((size_t)m * m * (n + 1), -1);
    auto solve = [&](auto&& self, int start, int end, int cnt) -> int{
        if (start > end) return 0;
        int& memo = dp[((size_t)start * m + end) * (n + 1) + cnt];
        if (memo != -1) return memo;

        int k = cnt + len[start], res = k * k + self(self, start + 1, end, 0);
        for (int mid = start + 1; mid <= end; mid++){
            if (color[start] == color[mid]) res = max(res, self(self, start + 1, mid - 1, 0) + self(self, mid, end, k));
        }
        return memo = res;
    };

    return solve(solve, 0, m - 1, 0);
}

int main(){
    /// The first block must skip its nearest equal neighbour to reach 50, random inputs this long rarely need that
    assert(blocks_max_score(vector<int>{0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0}) == 50);

    for (long long it = 0; it < stress::scaled(300); it++){
        /// skipping a nearer equal block first pays off at n = 12
        int n = stress::rand_int(0, it % 4 ? 9 : 12), colors = stress::rand_int(1, 3);
        vector<int> v(n);
        for (auto& x : v) x = stress::rand_int(1, colors) * (it % 2 ? 1 : 1000000);  /// compared, never indexed
        map<vector<int>, int> memo;
        assert(blocks_max_score(v) == brute(v, memo));
    }

    /// Every row of up to 8 blocks over 3 colors, so each carry chain shape near the cap is hit
    for (int n = 0; n <= 8; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 3;
        for (int code = 0; code < total; code++){
            vector<int> v(n);
            int rest = code;
            for (auto& x : v) x = rest % 3, rest /= 3;
            map<vector<int>, int> memo;
            assert(blocks_max_score(v) == brute(v, memo));
        }
    }

    /// Long runs of few colors make the carry reach max_carry exactly, where an off-by-one reads the next row
    for (long long it = 0; it < stress::scaled(60); it++){
        int n = stress::rand_int(1, it % 10 ? 60 : 200), colors = stress::rand_int(1, 4);
        int max_run = stress::rand_int(1, 8);
        vector<int> v;
        while ((int)v.size() < n){
            int color = stress::rand_int(1, colors), run = stress::rand_int(1, max_run);
            while (run-- && (int)v.size() < n) v.push_back(color);
        }
        assert(blocks_max_score(v) == memo_reference(v));
    }

    /// One color clears in a single click, distinct colors only one at a time, up to UVA's 200 block limit
    for (long long it = 0; it <= stress::scaled(10); it++){
        int n = it ? stress::rand_int(1, 60) : 200;
        vector<int> same(n, 7), distinct(n);
        iota(distinct.begin(), distinct.end(), 0);
        assert(blocks_max_score(same) == n * n);
        assert(blocks_max_score(distinct) == n);
    }

    /// 200 runs of 2 to 4 colors are the slowest rows, only bounded: >= 1 point per block, <= 1 click per color
    for (long long it = 0; it <= stress::scaled(2); it++){
        int n = 200, colors = it ? stress::rand_int(2, 4) : 2;
        vector<int> v(n);
        map<int, int> count;
        for (int i = 0; i < n; i++) v[i] = it ? stress::rand_int(1, colors) : i % 2, count[v[i]]++;

        int upper = 0;
        for (auto [color, c] : count) upper += c * c;
        int score = blocks_max_score(v);
        assert(n <= score && score <= upper);
    }

    return 0;
}
