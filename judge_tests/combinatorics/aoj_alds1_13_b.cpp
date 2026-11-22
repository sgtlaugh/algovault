// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_13_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/permutation_index.cpp"
#undef main

int main(){
    array<int, 9> start;
    for (int& x : start){
        if (scanf("%d", &x) != 1) return 0;
    }

    array<int, 9> goal = {1, 2, 3, 4, 5, 6, 7, 8, 0};
    PermutationIndex idx(9);
    vector<int> dist(362880, -1);
    queue<array<int, 9>> q;
    dist[idx.index(start)] = 0;
    q.push(start);

    int dx[4] = {-1, 1, 0, 0}, dy[4] = {0, 0, -1, 1};
    while (!q.empty()){
        array<int, 9> cur = q.front();
        q.pop();
        int d = dist[idx.index(cur)];
        if (cur == goal){
            printf("%d\n", d);
            return 0;
        }

        int blank = find(cur.begin(), cur.end(), 0) - cur.begin();
        for (int k = 0; k < 4; k++){
            int r = blank / 3 + dx[k], c = blank % 3 + dy[k];
            if (r < 0 || r >= 3 || c < 0 || c >= 3) continue;

            array<int, 9> next = cur;
            swap(next[blank], next[r * 3 + c]);
            long long id = idx.index(next);
            if (dist[id] == -1) dist[id] = d + 1, q.push(next);
        }
    }

    return 0;
}
