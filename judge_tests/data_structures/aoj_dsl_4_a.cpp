// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_4_A
// competitive-verifier: TLE 6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/interval_set.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<array<long long, 4>> rects(n);
    vector<long long> xs;
    for (auto& [x1, y1, x2, y2]: rects){
        if (scanf("%lld %lld %lld %lld", &x1, &y1, &x2, &y2) != 4) return 0;
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);
        xs.push_back(x1), xs.push_back(x2);
    }

    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());

    long long area = 0;
    for (int i = 0; i + 1 < (int)xs.size(); i++){
        IntervalSet<long long> slab;
        for (auto& [x1, y1, x2, y2]: rects){
            if (x1 <= xs[i] && xs[i + 1] <= x2) slab.add(y1, y2);
        }
        area += slab.covered_length() * (xs[i + 1] - xs[i]);
    }

    printf("%lld\n", area);
    return 0;
}
