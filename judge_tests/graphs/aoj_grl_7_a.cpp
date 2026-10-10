// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_7_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/matroid_intersection.cpp"
#undef main

int main(){
    int x, y, m;
    if (scanf("%d %d %d", &x, &y, &m) != 3) return 0;
    vector<int> left(m), right(m);
    for (int i = 0; i < m; i++){
        if (scanf("%d %d", &left[i], &right[i]) != 2) return 0;
    }

    /// An edge set is a matching iff no two edges share a left endpoint and no two share a right endpoint
    PartitionMatroid by_left(left, vector<int>(x, 1)), by_right(right, vector<int>(y, 1));
    printf("%d\n", (int)matroid_intersection(m, by_left, by_right).size());
    return 0;
}
