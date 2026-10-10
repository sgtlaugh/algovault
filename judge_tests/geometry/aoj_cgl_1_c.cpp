// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_1_C
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    Point p0, p1;
    int q;
    if (scanf("%lld %lld %lld %lld %d", &p0.x, &p0.y, &p1.x, &p1.y, &q) != 5) return 0;

    while (q--){
        Point p2;
        if (scanf("%lld %lld", &p2.x, &p2.y) != 2) return 0;
        long long c = cross(p0, p1, p2);
        if (c > 0) puts("COUNTER_CLOCKWISE");
        else if (c < 0) puts("CLOCKWISE");
        else if (on_segment(p2, p0, p1)) puts("ON_SEGMENT");
        else if (dot(p0, p1, p2) < 0) puts("ONLINE_BACK");
        else puts("ONLINE_FRONT");
    }
    return 0;
}
