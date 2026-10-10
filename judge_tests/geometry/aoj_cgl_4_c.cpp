// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_4_C
// competitive-verifier: ERROR 1e-5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<PointF> poly(n);
    for (auto& p : poly){
        if (scanf("%lf %lf", &p.x, &p.y) != 2) return 0;
    }

    int q;
    if (scanf("%d", &q) != 1) return 0;
    while (q--){
        PointF a, b;
        if (scanf("%lf %lf %lf %lf", &a.x, &a.y, &b.x, &b.y) != 4) return 0;
        vector<PointF> part = cut_polygon(poly, a, b);

        double area = 0;
        for (int i = 0, m = part.size(); i < m; i++){
            const PointF& p = part[i];
            const PointF& r = part[(i + 1) % m];
            area += p.x * r.y - p.y * r.x;
        }
        printf("%.10f\n", fabs(area) / 2);
    }
    return 0;
}
