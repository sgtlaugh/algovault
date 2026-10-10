// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0090
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    int n;
    while (scanf("%d", &n) == 1 && n != 0){
        vector<Point> centers(n);
        for (auto& p : centers){
            if (scanf("%Lf,%Lf", &p.x, &p.y) != 2) return 0;
        }

        /// A point lies under k unit stickers exactly when a unit circle around it covers k centers
        printf("%d\n", max_circle_cover(centers, 1).first);
    }
    return 0;
}
