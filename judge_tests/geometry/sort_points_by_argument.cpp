// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sort_points_by_argument
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
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
    int n = read_int();
    vector<Point> points(n);
    for (auto& p : points){
        p.x = read_int();
        p.y = read_int();
    }

    /// angle_less orders [0, 2 pi) with (0, 0) first, atan2 wants (-pi, pi]: lower half plane, then [0, pi), then pi
    auto group = [](const Point& p){ return p.y < 0 ? 0 : (p.y == 0 && p.x < 0 ? 2 : 1); };
    sort(points.begin(), points.end(), [&](const Point& a, const Point& b){
        int ga = group(a), gb = group(b);
        if (ga != gb) return ga < gb;
        return angle_less(a, b);
    });

    for (const Point& p : points) printf("%lld %lld\n", p.x, p.y);
    return 0;
}
