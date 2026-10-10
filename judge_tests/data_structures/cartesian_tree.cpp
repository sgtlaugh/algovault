// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/cartesian_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/monotonic_stack.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int();
    vector<int> a(n);
    for (auto& x: a) x = read_int();

    auto left = previous_smaller(a), right = next_smaller(a);

    /// The parent is the larger of the two nearest smaller neighbours, the global minimum has neither
    for (int i = 0; i < n; i++){
        int p = i;
        if (left[i] >= 0 && right[i] < n) p = a[left[i]] > a[right[i]] ? left[i] : right[i];
        else if (left[i] >= 0) p = left[i];
        else if (right[i] < n) p = right[i];

        printf("%d", p);
        putchar(i + 1 < n ? ' ' : '\n');
    }
    return 0;
}
