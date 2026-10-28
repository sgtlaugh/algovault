// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/nim_product_64
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/nimber.cpp"
#undef main
#define main fast_io_main
#include "../../code_library/misc/fast_io.cpp"
#undef main

int main(){
    int t;
    if (!fio::read(t)) return 0;

    while (t--){
        unsigned long long a, b;
        fio::read(a, b);
        fio::write(nim_product(a, b));
    }

    return 0;
}
