/***
 *
 * UVA 10559 Blocks - https://onlinejudge.org/external/105/10559.pdf
 * Maximum score for clearing a row of colored blocks, one click removing a maximal run of k equal blocks for k^2 points
 *
 * Complexity: O(m^3 n) time, O(m^2 n) memory for n blocks in m maximal runs (32 MB at n = 200)
 *
 * blocks_max_score(blocks): the maximum score, 0 for an empty row, blocks any random-access container compared with ==
 *
 * Notes:
 *   - A similar problem solved the same way - https://www.spoj.com/problems/DINGRP/
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
int blocks_max_score(const Sequence& blocks){
    int n = blocks.size();
    vector<int> run_start, run_len;

    for (int i = 0; i < n; i++){
        if (!run_len.empty() && blocks[run_start.back()] == blocks[i]) run_len.back()++;
        else run_start.push_back(i), run_len.push_back(1);
    }

    /// cnt + blocks in runs start..end <= n and the range is non-empty, so cnt < n
    int m = run_len.size();
    vector<int> dp((size_t)m * m * n, -1);

    auto solve = [&](auto&& self, int start, int end, int cnt) -> int{
        if (start > end) return 0;
        int& memo = dp[((size_t)start * m + end) * n + cnt];
        if (memo != -1) return memo;

        int k = cnt + run_len[start], res = k * k + self(self, start + 1, end, 0);
        for (int mid = start + 1; mid <= end; mid++){
            if (blocks[run_start[start]] == blocks[run_start[mid]]){
                res = max(res, self(self, start + 1, mid - 1, 0) + self(self, mid, end, k));
            }
        }

        return memo = res;
    };

    return solve(solve, 0, m - 1, 0);
}

int main(){
    assert(blocks_max_score(vector<int>{1, 2, 2, 2, 2, 3, 3, 3, 1}) == 29);
    assert(blocks_max_score(vector<int>{1}) == 1);
    assert(blocks_max_score(vector<int>{1, 2, 2, 2, 2, 2, 3, 2}) == 38);
    assert(blocks_max_score(vector<int>{}) == 0);
    assert(blocks_max_score(vector<int>{4, 4, 4}) == 9);
    assert(blocks_max_score(vector<int>{1, 2, 3}) == 3);
    assert(blocks_max_score(string("aabaa")) == 17);

    return 0;
}
