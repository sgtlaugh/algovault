// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_I
// competitive-verifier: ERROR 1e-6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Circle c[2];
    for (auto& ci : c){
        long long x, y, r;
        if (scanf("%lld %lld %lld", &x, &y, &r) != 3) return 0;
        ci = {{(Float)x, (Float)y}, (Float)r};
    }

    printf("%.12Lf\n", circle_circle_area(c[0], c[1]));
    return 0;
}
