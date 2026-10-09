#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/xor_segment_tree.cpp"
#undef main

long long brute_query(const vector<long long>& a, int l, int r, int x){
    long long sum = 0;
    for (int p = l; p <= r; p++) sum += a[p ^ x];
    return sum;
}

/// compares against relabelling every index p -> p ^ x and summing directly, after random point adds
int main(){
    for (int k = 0; k <= 3; k++){
        int n = 1 << k;
        vector<long long> a(n);
        for (auto& v : a) v = stress::rand_int(-100, 100);

        XorSegmentTree<long long> st(a);
        for (int round = 0; round < 4; round++){
            for (int x = 0; x < n; x++){
                for (int l = 0; l < n; l++){
                    for (int r = l; r < n; r++) assert(st.query(l, r, x) == brute_query(a, l, r, x));
                }
            }

            int i = stress::rand_int(0, n - 1);
            long long v = stress::rand_int(-100, 100);
            st.add(i, v);
            a[i] += v;
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = 1 << stress::rand_int(0, it % 50 ? 9 : 14);
        long long lim = it % 3 ? 1000 : 1000000000000LL;
        vector<long long> a(n);
        for (auto& v : a) v = stress::rand_int(-lim, lim);

        XorSegmentTree<long long> st(a);
        for (int q = 0; q < 60; q++){
            if (q % 3 == 0){
                int i = stress::rand_int(0, n - 1);
                long long v = stress::rand_int(-lim, lim);
                st.add(i, v);
                a[i] += v;
                continue;
            }

            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1), x = stress::rand_int(0, n - 1);
            if (q % 7 == 0) l = 0, r = n - 1;
            assert(st.query(l, r, x) == brute_query(a, l, r, x));
        }
    }

    return 0;
}
