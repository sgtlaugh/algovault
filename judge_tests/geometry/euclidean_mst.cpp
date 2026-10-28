// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/euclidean_mst
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/delaunay.cpp"
#undef main

#define main mst_main
#include "../../code_library/graphs/minimum_spanning_tree.cpp"
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

    Delaunay dt(points);
    Kruskal mst(n);
    for (auto [u, v] : dt.edges){
        long long dx = points[u].x - points[v].x, dy = points[u].y - points[v].y;
        mst.add_edge(u, v, dx * dx + dy * dy);
    }

    for (int i : mst.solve().chosen) printf("%d %d\n", dt.edges[i][0], dt.edges[i][1]);
    return 0;
}
