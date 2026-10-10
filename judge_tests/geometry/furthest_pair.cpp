// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/furthest_pair
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/convex_hull.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != '-' && (c < '0' || c > '9')) c = getchar();

    bool neg = c == '-';
    if (neg) c = getchar();
    long long v = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) v = v * 10 + c - '0';
    return neg ? -v : v;
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

        array<Point, 2> best = farthest_pair(get_convex_hull(points));
        int i = find(points.begin(), points.end(), best[0]) - points.begin();
        int j = find(points.begin(), points.end(), best[1]) - points.begin();
        if (i == j) i = 0, j = 1;
        printf("%d %d\n", i, j);
    }
    return 0;
}
