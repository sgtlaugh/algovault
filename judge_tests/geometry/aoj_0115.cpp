// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0115
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry_3D.cpp"
#undef main

int main(){
    Point3 p[5];
    for (auto& q : p){
        if (scanf("%lf %lf %lf", &q.x, &q.y, &q.z) != 3) return 0;
    }

    /// Integer coordinates up to 100: a beam that misses passes far more than 1e-9 from the barrier
    bool blocked = dist_segment_triangle(p[0], p[1], p[2], p[3], p[4]) < 1e-9;
    puts(blocked ? "MISS" : "HIT");
    return 0;
}
