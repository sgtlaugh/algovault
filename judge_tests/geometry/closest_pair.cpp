// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/closest_pair
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/closest_pair.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != '-' && (c < '0' || c > '9')) c = getchar();
    bool neg = c == '-';
    if (neg) c = getchar();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + (c - '0');
    return neg ? -x : x;
}

int main(){
    int t = read_int();
    while (t--){
        int n = read_int();
        vector<Point> points(n);
        for (auto& p : points){
            p.x = read_int();
            p.y = read_int();
        }

        ClosestPairResult res = closest_pair(points);
        printf("%d %d\n", res.i, res.j);
    }
    return 0;
}
