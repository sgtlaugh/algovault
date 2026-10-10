// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_5_A
// competitive-verifier: ERROR 1e-6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/closest_pair.cpp"
#undef main

/// Reads a decimal with at most 6 fractional digits exactly as an integer scaled by 1e6
long long read_scaled(){
    int c = getchar();
    while (c != '-' && c != '.' && (c < '0' || c > '9')) c = getchar();

    bool neg = c == '-';
    if (neg) c = getchar();
    long long v = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) v = v * 10 + c - '0';

    int digits = 0;
    if (c == '.'){
        for (c = getchar(); c >= '0' && c <= '9' && digits < 6; c = getchar(), digits++) v = v * 10 + c - '0';
    }
    for (; digits < 6; digits++) v *= 10;
    return neg ? -v : v;
}

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<Point> points(n);
    for (auto& p : points){
        p.x = read_scaled();
        p.y = read_scaled();
    }

    printf("%.10Lf\n", sqrtl((long double)closest_pair(points).dist2) / 1e6L);
    return 0;
}
