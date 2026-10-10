// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/many_aplusb
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/fast_io.cpp"
#undef main

int main(){
    int t;
    if (!fio::read(t)) return 0;

    while (t--){
        long long a, b;
        fio::read(a, b);
        fio::write(a + b);
    }
    return 0;
}
