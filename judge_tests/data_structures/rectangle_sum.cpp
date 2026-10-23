// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/rectangle_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

int main(){
    int n, q;
    scanf("%d %d", &n, &q);

    vector<array<int, 3>> points(n);
    vector<int> ys(n);
    for (int i = 0; i < n; i++){
        scanf("%d %d %d", &points[i][0], &points[i][1], &points[i][2]);
        ys[i] = points[i][1];
    }
    sort(points.begin(), points.end());
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    /// Each query splits into prefix-in-x events: +rect(x < r) - rect(x < l)
    vector<array<int, 5>> events;
    for (int i = 0; i < q; i++){
        int l, d, r, u;
        scanf("%d %d %d %d", &l, &d, &r, &u);
        int lo = lower_bound(ys.begin(), ys.end(), d) - ys.begin();
        int hi = lower_bound(ys.begin(), ys.end(), u) - ys.begin() - 1;
        events.push_back({r, i, lo, hi, 1});
        events.push_back({l, i, lo, hi, -1});
    }
    sort(events.begin(), events.end());

    FenwickPointUpdate<long long> fen(ys.size());
    vector<long long> res(q, 0);
    size_t j = 0;
    for (auto& e: events){
        for (; j < points.size() && points[j][0] < e[0]; j++){
            int y = lower_bound(ys.begin(), ys.end(), points[j][1]) - ys.begin();
            fen.update(y, points[j][2]);
        }
        res[e[1]] += e[4] * fen.query(e[2], e[3]);
    }

    for (auto x: res) printf("%lld\n", x);
    return 0;
}
