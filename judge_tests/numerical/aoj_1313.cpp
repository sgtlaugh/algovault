// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/1313
// competitive-verifier: TLE 8
// competitive-verifier: ERROR 1e-3
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/numerical/numerical_integration.cpp"
#undef main

/// Length of the convex polygon's cut by the vertical line at x, 0 when the line misses it
double width(const vector<pair<int, int>>& poly, double x){
    double lo = 1e18, hi = -1e18;
    for (size_t i = 0; i < poly.size(); i++){
        auto [x1, y1] = poly[i];
        auto [x2, y2] = poly[(i + 1) % poly.size()];
        if (x < min(x1, x2) || x > max(x1, x2)) continue;

        if (x1 == x2) lo = min({lo, (double)y1, (double)y2}), hi = max({hi, (double)y1, (double)y2});
        else{
            double y = y1 + (double)(y2 - y1) * (x - x1) / (x2 - x1);
            lo = min(lo, y), hi = max(hi, y);
        }
    }

    return hi > lo ? hi - lo : 0;
}

int main(){
    int m, n;
    while (scanf("%d %d", &m, &n) == 2 && (m || n)){
        vector<pair<int, int>> a(m), b(n);
        int lo_a = INT_MAX, hi_a = INT_MIN, lo_b = INT_MAX, hi_b = INT_MIN;
        vector<int> xs;
        for (auto& [x, y] : a){
            if (scanf("%d %d", &x, &y) != 2) return 0;
            xs.push_back(x), lo_a = min(lo_a, x), hi_a = max(hi_a, x);
        }
        for (auto& [x, z] : b){
            if (scanf("%d %d", &x, &z) != 2) return 0;
            xs.push_back(x), lo_b = min(lo_b, x), hi_b = max(hi_b, x);
        }
        int lo = max(lo_a, lo_b), hi = min(hi_a, hi_b);
        xs.push_back(lo), xs.push_back(hi);
        sort(xs.begin(), xs.end());
        xs.erase(unique(xs.begin(), xs.end()), xs.end());

        /// Between consecutive vertex x inside both spans the cross section is a product of two linear widths
        double volume = 0;
        auto area = [&](double x){ return width(a, x) * width(b, x); };
        for (size_t i = 0; i + 1 < xs.size(); i++){
            if (lo <= xs[i] && xs[i + 1] <= hi) volume += simpson(xs[i], xs[i + 1], area, 1);
        }

        printf("%.10f\n", volume);
    }

    return 0;
}
