// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_2_D
// competitive-verifier: ERROR 1e-8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    while (q--){
        Point p[4];
        for (auto& pt : p){
            if (scanf("%lld %lld", &pt.x, &pt.y) != 2) return 0;
        }

        double res = 0;
        if (!segments_intersect(p[0], p[1], p[2], p[3])){
            res = min({dist_point_segment(p[0], p[2], p[3]), dist_point_segment(p[1], p[2], p[3]),
                       dist_point_segment(p[2], p[0], p[1]), dist_point_segment(p[3], p[0], p[1])});
        }
        printf("%.10f\n", res);
    }
    return 0;
}
