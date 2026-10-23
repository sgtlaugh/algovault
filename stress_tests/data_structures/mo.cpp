#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
#undef main

/// Distinct values and sum of every range, empty ones included, against direct scans; get also sees the right query index
/// check_calls holds add and remove to the header bound with constant 1, which random inputs stay under and a too small block exceeds
void check_array(int n, int q, int max_v, int verify = 1 << 30, bool check_calls = false){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(0, max_v);
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        l = stress::rand_int(0, n);  /// r = l - 1 is an empty range, l = n included
        r = l < n && stress::rand_int(0, 4) ? stress::rand_int(l, n - 1) : l - 1;
    }

    vector<int> freq(max_v + 1, 0);
    int distinct = 0, window = 0;
    long long sum = 0;

    for (int round = 0; round < 2; round++){
        long long calls = 0;
        auto got = mo(n, queries,
            [&](int i){ distinct += freq[a[i]]++ == 0, sum += a[i], window++, calls++; },
            [&](int i){ distinct -= --freq[a[i]] == 0, sum -= a[i], window--, calls++; },
            [&](int qi){ return make_tuple(distinct, sum, qi); });
        assert(window == 0 && sum == 0 && distinct == 0 && (int)got.size() == q);
        assert(!check_calls || calls <= n * sqrt((double)q) + q);
        for (int i = 0; i < q && i < verify; i++){
            auto [l, r] = queries[i];
            auto [got_distinct, got_sum, got_qi] = got[i];
            assert(got_qi == i);
            assert(got_sum == accumulate(a.begin() + l, a.begin() + r + 1, 0LL));
            assert(got_distinct == (int)set<int>(a.begin() + l, a.begin() + r + 1).size());
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_array(stress::rand_int(1, 60), stress::rand_int(1, 60), it % 2 ? 3 : 50);
    }

    check_array(100000, 100000, 1000, 300, true);
    check_array(100000, 1000, 1000, 300, true);

    return 0;
}
