// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_6_A
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/combsort.cpp"
#undef main
#define main fast_io_main
#include "../../code_library/misc/fast_io.cpp"
#undef main

int main(){
    vector<int> a;
    if (!fio::read(a)) return 0;

    combsort(a.begin(), a.end());
    fio::write(a);
    return 0;
}
