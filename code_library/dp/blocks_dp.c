/***
 *
 * UVA 10559 Blocks - https://onlinejudge.org/external/105/10559.pdf
 * Given an array of blocks, need to maximize the score following the rules below
 * Color of any block i is stored in blocks[i]
 * Blocks can be clicked on and clicking a block i clears all consecutive blocks to the left and right of i
 * If k blocks are cleared in one move, score is incremented by k * k
 * So if the blocks are 1, 2, 2, 2, 2, 2, 3, 2, then after clicking on the fourth box it becomes 1, 3, 2 and a score of 25 is added
 *
 * It's one kind of a variation of matrix chain multiplication and useful in problems where subsequences need to be chained
 * A similar problem which can be solved using this method - https://www.spoj.com/problems/DINGRP/
 *
***/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#define MAX 201

int dp[MAX][MAX][MAX], run_color[MAX], run_len[MAX];

/// Best score for runs start..end when cnt blocks of run_color[start] wait to be cleared together with run start
int solve(int start, int end, int cnt){
    if (start > end) return 0;
    if (dp[start][end][cnt] != -1) return dp[start][end][cnt];

    int mid, k = cnt + run_len[start], res = k * k + solve(start + 1, end, 0);
    for (mid = start + 1; mid <= end; mid++){
        if (run_color[start] == run_color[mid]){
            int cur = solve(start + 1, mid - 1, 0) + solve(mid, end, k);
            if (cur > res) res = cur;
        }
    }

    return dp[start][end][cnt] = res;
}

int get_max_score(const int n, const int blocks[]){
    int i, j, m = 0;
    for (i = 0; i < n; i++){
        if (m && run_color[m - 1] == blocks[i]) run_len[m - 1]++;
        else run_color[m] = blocks[i], run_len[m++] = 1;
    }

    /// Only cnt < n is reachable, clearing just that instead of all 32 MB keeps multi-test inputs fast
    for (i = 0; i < m; i++){
        for (j = 0; j < m; j++) memset(dp[i][j], -1, sizeof(int) * n);
    }

    return solve(0, m - 1, 0);
}

int main(){
    const int n = 9;
    const int blocks[] = {1, 2, 2, 2, 2, 3, 3, 3, 1};
    assert(get_max_score(n, blocks) == 29);
    return 0;
}
