// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/segment_add_get_min
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/li_chao_tree.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    const long long limit = 1e9;
    LiChaoTree<long long> tree(-limit, limit);
    for (int i = 0; i < n; i++){
        long long l, r, a, b;
        if (scanf("%lld %lld %lld %lld", &l, &r, &a, &b) != 4) return 0;
        tree.add_segment(a, b, l, r - 1);
    }

    while (q--){
        int type;
        if (scanf("%d", &type) != 1) return 0;
        if (type == 0){
            long long l, r, a, b;
            if (scanf("%lld %lld %lld %lld", &l, &r, &a, &b) != 4) return 0;
            tree.add_segment(a, b, l, r - 1);
        }
        else{
            long long p;
            if (scanf("%lld", &p) != 1) return 0;
            long long y = tree.query(p);
            if (y == tree.NONE) puts("INFINITY");
            else printf("%lld\n", y);
        }
    }
    return 0;
}
