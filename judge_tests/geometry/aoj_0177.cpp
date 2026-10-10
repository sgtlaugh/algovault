// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0177
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    long double a, b, c, d;
    while (scanf("%Lf %Lf %Lf %Lf", &a, &b, &c, &d) == 4){
        if (a == -1 && b == -1 && c == -1 && d == -1) break;
        printf("%lld\n", llroundl(great_circle_distance(a, b, c, d, 6378.1L)));
    }
    return 0;
}
