#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/mo_with_updates.cpp"
#undef main

/// Hashes (position, value) so a window sum of hashes pins down exactly which cells, at which version, are inside
unsigned long long cell_hash(int i, int v){
    unsigned long long x = (unsigned long long)i * 1000003 + v + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/// Distinct values and cell hash of every timed range, empty ones included, against replaying the updates on a fresh copy
/// check_calls holds the callback count to 2 * the header bound, random inputs reach about 1.16 at large n
void check_updates(int n, int q, int u, int max_v, int verify = 1 << 30, bool check_calls = false){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(0, max_v);
    vector<pair<int, int>> updates(u);
    for (auto& [p, v] : updates) p = stress::rand_int(0, n - 1), v = stress::rand_int(0, max_v);
    vector<array<int, 3>> queries(q);
    for (auto& [l, r, t] : queries){
        l = stress::rand_int(0, n);  /// r = l - 1 is an empty range, l = n included
        r = l < n && stress::rand_int(0, 4) ? stress::rand_int(l, n - 1) : l - 1;
        t = stress::rand_int(0, u);
    }
    const vector<int> initial = a;

    vector<int> freq(max_v + 1, 0);
    int distinct = 0;
    unsigned long long hash = 0;
    double block = cbrt((double)n * n);

    for (int round = 0; round < 2; round++){
        long long calls = 0;
        auto got = mo_with_updates(a, queries, updates,
            [&](int i){ distinct += freq[a[i]]++ == 0, hash += cell_hash(i, a[i]), calls++; },
            [&](int i){ distinct -= --freq[a[i]] == 0, hash -= cell_hash(i, a[i]), calls++; },
            [&](int qi){ return make_tuple(distinct, hash, qi); });
        assert(distinct == 0 && hash == 0 && a == initial && (int)got.size() == q);
        assert(!check_calls || calls <= 2 * ((q + u) * block + block * block));

        for (int i = 0; i < q && i < verify; i++){
            auto [l, r, t] = queries[i];
            vector<int> b = initial;
            for (int j = 0; j < t; j++) b[updates[j].first] = updates[j].second;

            unsigned long long h = 0;
            for (int k = l; k <= r; k++) h += cell_hash(k, b[k]);
            auto [got_distinct, got_hash, got_qi] = got[i];
            assert(got_qi == i && got_hash == h && got_distinct == (int)set<int>(b.begin() + l, b.begin() + r + 1).size());
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_updates(stress::rand_int(1, 60), stress::rand_int(1, 60), stress::rand_int(0, 30), it % 2 ? 3 : 50);
    }

    check_updates(20000, 20000, 20000, 1000, 200, true);
    return 0;
}
