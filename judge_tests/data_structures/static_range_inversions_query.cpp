// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_inversions_query
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
#undef main

#define main fenwick_tree_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    vector<int> a(n);
    for (int& x: a){
        if (scanf("%d", &x) != 1) return 0;
    }

    vector<int> values = a;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    int m = values.size();
    for (int& x: a) x = lower_bound(values.begin(), values.end(), x) - values.begin() + 1;

    vector<pair<int, int>> queries(q);
    for (auto& [l, r]: queries){
        if (scanf("%d %d", &l, &r) != 2) return 0;
        r--;
    }

    /// The callbacks only get a position, so the window [left, right] tells which end moved
    FenwickPointUpdate<int> counts(m);
    int left = 0, right = -1;
    long long inversions = 0;
    vector<long long> answer(q);

    auto add = [&](int i){
        int size = right - left + 1;
        if (i == left - 1) inversions += counts.query(a[i] - 1), left--;
        else inversions += size - counts.query(a[i]), right++;
        counts.update(a[i], 1);
    };
    auto remove = [&](int i){
        counts.update(a[i], -1);
        int size = right - left;
        if (i == left) inversions -= counts.query(a[i] - 1), left++;
        else inversions -= size - counts.query(a[i]), right--;
    };

    mo(n, queries, add, remove, [&](int qi){ answer[qi] = inversions; });

    for (long long x: answer) printf("%lld\n", x);
    return 0;
}
