// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/area_of_union_of_rectangles
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/rectangle_union.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + (c - '0');
    return x;
}

int main(){
    int n = read_int();
    vector<Rectangle> rects(n);
    for (auto& r : rects){
        r.x1 = read_int();
        r.y1 = read_int();
        r.x2 = read_int();
        r.y2 = read_int();
    }

    printf("%lld\n", rectangle_union_area(rects));
    return 0;
}
