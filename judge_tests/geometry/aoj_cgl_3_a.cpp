// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_3_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<Point> poly(n);
    for (auto& p : poly){
        if (scanf("%lld %lld", &p.x, &p.y) != 2) return 0;
    }

    long long a = llabs(area2(poly));
    printf("%lld.%d\n", a / 2, a % 2 ? 5 : 0);
    return 0;
}
