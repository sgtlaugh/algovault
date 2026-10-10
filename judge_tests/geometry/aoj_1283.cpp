// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/1283
// competitive-verifier: ERROR 1e-5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/halfplane_intersection.cpp"
#undef main

int main(){
    int n;
    while (scanf("%d", &n) == 1 && n != 0){
        vector<Point> poly(n);
        for (auto& p : poly){
            if (scanf("%Lf %Lf", &p.x, &p.y) != 2) return 0;
        }

        /// A point at distance d from the sea survives every edge pushed inward by d along its left normal
        auto fits = [&](long double d){
            vector<Halfplane> h;
            for (int i = 0; i < n; i++){
                Point a = poly[i], b = poly[(i + 1) % n], dir = b - a;
                Point shift = Point(-dir.y, dir.x) * (d / hypotl(dir.x, dir.y));
                h.emplace_back(a + shift, b + shift);
            }
            return !halfplane_intersection(h, 1e5).empty();
        };

        long double lo = 0, hi = 1e4;
        for (int it = 0; it < 60; it++){
            long double mid = (lo + hi) / 2;
            if (fits(mid)) lo = mid;
            else hi = mid;
        }
        printf("%.6Lf\n", lo);
    }
    return 0;
}
