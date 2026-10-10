// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/queue_operate_all_composite
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/sliding_window_aggregation.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    long long a, b;
};

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    auto op = [](const Affine& f, const Affine& g){ return Affine{f.a * g.a % MOD, (f.b * g.a + g.b) % MOD}; };
    auto swag = SlidingWindowAggregation<Affine, decltype(op)>(Affine{1, 0}, op);
    while (q--){
        int t;
        if (scanf("%d", &t) != 1) return 0;
        if (t == 0){
            long long a, b;
            if (scanf("%lld %lld", &a, &b) != 2) return 0;
            swag.push({a, b});
        }
        else if (t == 1) swag.pop();
        else{
            long long x;
            if (scanf("%lld", &x) != 1) return 0;
            Affine f = swag.fold();
            printf("%lld\n", (f.a * x + f.b) % MOD);
        }
    }
    return 0;
}
