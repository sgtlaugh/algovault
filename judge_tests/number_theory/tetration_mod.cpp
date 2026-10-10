// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/tetration_mod
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/power_tower.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long a, b, m;
        if (scanf("%lld %lld %lld", &a, &b, &m) != 3) return 0;

        /// 0^0 = 1 makes the towers of zeros alternate 1, 0, 1, 0, ... and the library requires every a[i] >= 1
        if (b == 0) printf("%lld\n", 1 % m);
        else if (a == 0) printf("%lld\n", (b % 2 == 0) % m);
        else{
            /// The phi chain of m <= 1e9 reaches 1 within about 60 levels, so copies past 128 never change the result
            vector<long long> tower(min(b, 128LL), a);
            printf("%lld\n", power_tower(tower, m));
        }
    }
    return 0;
}
