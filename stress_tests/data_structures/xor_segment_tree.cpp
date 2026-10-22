#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/xor_segment_tree.cpp"
#undef main

struct Node{
    long long sum, pre, suf, best;
};

struct Join{
    Node operator()(const Node& a, const Node& b) const{
        return {a.sum + b.sum, max(a.pre, a.sum + b.pre), max(b.suf, b.sum + a.suf), max({a.best, b.best, a.suf + b.pre})};
    }
};

long long brute_query(const vector<long long>& a, int l, int r, int x){
    long long sum = 0;
    for (int p = l; p <= r; p++) sum += a[p ^ x];
    return sum;
}

/// Kadane over the relabelled sequence, with the empty subarray allowed
long long brute_best(const vector<long long>& a, int l, int r, int x){
    long long best = 0, cur = 0;
    for (int p = l; p <= r; p++) cur = max(0LL, cur + a[p ^ x]), best = max(best, cur);
    return best;
}

/// Every query against relabelling each index p -> p ^ x directly
int main(){
    /// Concatenation and max subarray sum do not commute, so a fold out of p order shows
    for (int k = 0; k <= 4; k++){
        int n = 1 << k;
        vector<string> words(n);
        for (int i = 0; i < n; i++) words[i] = string(1, 'a' + i);

        XorSegmentTree concat(words, [](const string& x, const string& y){ return x + y; });
        for (int x = 0; x < n; x++){
            for (int l = 0; l < n; l++){
                string expected;
                for (int r = l; r < n; r++){
                    expected += words[r ^ x];
                    assert(concat.query(l, r, x) == expected);
                }
            }
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = 1 << stress::rand_int(0, it % 50 ? 7 : 12);
        vector<long long> a(n);
        vector<Node> nodes(n);
        for (int i = 0; i < n; i++){
            a[i] = stress::rand_int(-10, 10);
            long long m = max(0LL, a[i]);
            nodes[i] = {a[i], m, m, m};
        }

        XorSegmentTree<Node, Join> st(nodes);
        for (int q = 0; q < 30; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1), x = stress::rand_int(0, n - 1);
            assert(st.query(l, r, x).best == brute_best(a, l, r, x));
        }
    }

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
