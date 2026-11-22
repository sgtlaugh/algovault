// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0125
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/dates.cpp"
#undef main

int main(){
    long long y1, m1, d1, y2, m2, d2;
    while (scanf("%lld %lld %lld %lld %lld %lld", &y1, &m1, &d1, &y2, &m2, &d2) == 6){
        if (y1 < 0 || m1 < 0 || d1 < 0 || y2 < 0 || m2 < 0 || d2 < 0) break;
        printf("%lld\n", day_number(y2, m2, d2) - day_number(y1, m1, d1));
    }

    return 0;
}
