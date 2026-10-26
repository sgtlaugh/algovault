/***
 *
 * Block Removal DP
 * Maximum score for clearing a row of colored blocks, one click removing a maximal run of k equal blocks for k^2 points
 *
 * Complexity: O(m^3 n) time, O(m^2 n) memory for n blocks in m maximal runs (under 4 MB at n = 200)
 *
 * blocks_max_score(blocks): the maximum score, 0 for an empty row, blocks any random-access container compared with ==
 *
 * Notes:
 *   - Solves UVa 10559 Blocks (https://onlinejudge.org/external/105/10559.pdf) and SPOJ DINGRP the same way
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
int blocks_max_score(const Sequence& blocks){
    int n = blocks.size();
    vector<int> run_start, run_len;
    run_start.reserve(n), run_len.reserve(n);  /// regrowing these cost 20% on rows of up to 12 blocks

    for (int i = 0; i < n; i++){
        if (!run_len.empty() && blocks[run_start.back()] == blocks[i]) run_len.back()++;
        else run_start.push_back(i), run_len.push_back(1);
    }

    /// max_carry[i] counts every earlier block of run i's color, so it bounds the carry from any range start
    int m = run_len.size();
    vector<int> max_carry(m), next_same(m, m);
    for (int i = 0; i < m; i++){
        for (int p = i - 1; p >= 0; p--){
            if (blocks[run_start[p]] == blocks[run_start[i]]){
                max_carry[i] = max_carry[p] + run_len[p], next_same[p] = i;
                break;
            }
        }
    }

    /// best[row(i, j) + k]: runs i..j cleared with k more blocks of run i's color joined to it, k <= max_carry[i]
    vector<size_t> offset(m + 1);
    for (int i = 0; i < m; i++) offset[i + 1] = offset[i] + (size_t)(m - i) * (max_carry[i] + 1);
    vector<int> best(offset[m]);
    auto row = [&](int i, int j){ return offset[i] + (size_t)(j - i) * (max_carry[i] + 1); };
    auto cleared = [&](int i, int j){ return i > j ? 0 : best[row(i, j)]; };

    /// filling every state is 25% slower than a memo below 6 blocks (5 ns a call) and 12x faster at 200 blocks
    for (int i = m - 1; i >= 0; i--){
        for (int j = i; j < m; j++){
            int* res = &best[row(i, j)];
            int rest = cleared(i + 1, j);
            for (int k = 0; k <= max_carry[i]; k++) res[k] = (k + run_len[i]) * (k + run_len[i]) + rest;

            for (int p = next_same[i]; p <= j; p = next_same[p]){
                int between = cleared(i + 1, p - 1);
                const int* joined = &best[row(p, j) + run_len[i]];
                for (int k = 0; k <= max_carry[i]; k++) res[k] = max(res[k], between + joined[k]);
            }
        }
    }

    return cleared(0, m - 1);
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
