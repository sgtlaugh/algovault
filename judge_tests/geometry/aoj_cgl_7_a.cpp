// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Circle a, b;
    if (scanf("%Lf %Lf %Lf %Lf %Lf %Lf", &a.center.x, &a.center.y, &a.r, &b.center.x, &b.center.y, &b.r) != 6) return 0;

    printf("%d\n", (int)common_tangents(a, b).size());
    return 0;
}
