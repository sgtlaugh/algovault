// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/two_square_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/sum_of_two_squares.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    TwoSquares squares;
    while (q--){
        long long n;
        if (scanf("%lld", &n) != 1) return 0;
        auto pairs = squares.pairs(n);
        printf("%d\n", (int)pairs.size());
        for (auto [x, y] : pairs) printf("%lld %lld\n", x, y);
    }
    return 0;
}
