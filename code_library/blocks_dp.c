#include <stdio.h>
#include <string.h>

#define MAX 201

int dp[MAX][MAX][MAX];

int solve(const int n, const int blocks[], int start, int end, int cnt){
    if (start > end) return 0;
    if (dp[start][end][cnt] != -1) return dp[start][end][cnt];

    int mid;
    int res = (cnt + 1) * (cnt + 1) + solve(n, blocks, start + 1, end, 0);

    for (mid = start + 1; mid <= end; mid++){
        if (blocks[start] == blocks[mid]){
            int x = solve(n, blocks, start + 1, mid - 1, 0) + solve(n, blocks, mid, end, cnt + 1);
            if (x > res) res = x;
        }
    }

    return (dp[start][end][cnt] = res);
}

/***
 *
 * UVA 10559 Blocks - https://onlinejudge.org/external/105/10559.pdf
 * given an array of blocks, need to maximize the score following the rules below
 * color of any block i is stored in blocks[i]
 * boxes can be clicked on and clicking a box i clears all consecutive boxes to the left and right of i
 * if k blocks are cleared in one move, score is incremented by k * k
 * so if the blocks are 1, 2, 2, 2, 2, 2, 3, 2, then after clicking on the fourth box it becomes 1, 3, 2 and a score of 25 is added
 *
 * it's one kind of a variation of matrix chain multiplication and useful in problems where subsequences need to be chained
 * a similar problem which can be solved using this method - https://www.spoj.com/problems/DINGRP/
 *
***/

int get_max_score(const int n, const int blocks[]){
    memset(dp, -1, sizeof(dp));
    return solve(n, blocks, 0, n - 1, 0);
}

int main(){
    const int n = 9;
    const int blocks[] = {1, 2, 2, 2, 2, 3, 3, 3, 1};

    printf("%d\n", get_max_score(n, blocks));
    return 0;
}
