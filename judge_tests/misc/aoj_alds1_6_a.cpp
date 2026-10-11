// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_6_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/combsort.cpp"
#undef main
#define main fast_io_main
#include "../../code_library/misc/fast_io.cpp"
#undef main

int main(){
    int n;
    if (!fio::read(n)) return 0;
    vector<int> a(n);
    fio::read(a);

    combsort(a.begin(), a.end());
    fio::write(a);
    return 0;
}
