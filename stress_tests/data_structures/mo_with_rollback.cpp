#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/mo_with_rollback.cpp"
#undef main

/// Hashes (position, value) so a window sum of hashes pins down exactly which cells, at which version, are inside
unsigned long long cell_hash(int i, int v){
    unsigned long long x = (unsigned long long)i * 1000003 + v + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/// max(value * count) and cell hash of every range, empty ones included, against direct scans, with add-only state and a contiguity check
/// check_calls holds the add count to 2 * the header bound, random inputs reach about 0.99 at large n
void check_rollback(int n, int q, int max_v, int verify = 1 << 30, bool check_calls = false){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(1, max_v);
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        l = stress::rand_int(0, n);  /// r = l - 1 is an empty range, l = n included
        r = l < n && stress::rand_int(0, 4) ? stress::rand_int(l, n - 1) : l - 1;
    }

    struct Saved{
        int pos, wl, wr;
        long long best;
        unsigned long long hash;
    };
    vector<int> freq(max_v + 1, 0);
    vector<Saved> history;
    vector<size_t> marks;
    int wl = 0, wr = -1;
    long long best = 0;
    unsigned long long hash = 0;

    for (int round = 0; round < 2; round++){
        long long calls = 0;
        auto got = mo_with_rollback(n, queries,
            [&](int i){
                assert(wl > wr || i == wl - 1 || i == wr + 1);
                calls++;
                history.push_back({i, wl, wr, best, hash});
                if (wl > wr) wl = wr = i;
                else wl = min(wl, i), wr = max(wr, i);
                best = max(best, (long long)a[i] * ++freq[a[i]]), hash += cell_hash(i, a[i]);
            },
            [&](){
                assert(marks.size() < 2);
                marks.push_back(history.size());
            },
            [&](){
                assert(!marks.empty());
                for (; history.size() > marks.back(); history.pop_back()){
                    Saved s = history.back();
                    freq[a[s.pos]]--, wl = s.wl, wr = s.wr, best = s.best, hash = s.hash;
                }
                marks.pop_back();
            },
            [&](int qi){ return make_tuple(best, hash, qi); });
        assert(marks.empty() && history.empty() && wl > wr && best == 0 && hash == 0 && (int)got.size() == q);
        assert(!check_calls || calls <= 2 * (n * sqrt((double)q) + q));

        for (int i = 0; i < q && i < verify; i++){
            auto [l, r] = queries[i];
            vector<int> count(max_v + 1, 0);
            long long expected = 0;
            unsigned long long h = 0;
            for (int k = l; k <= r; k++) expected = max(expected, (long long)a[k] * ++count[a[k]]), h += cell_hash(k, a[k]);
            auto [got_best, got_hash, got_qi] = got[i];
            assert(got_qi == i && got_best == expected && got_hash == h);
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_rollback(stress::rand_int(1, 60), stress::rand_int(1, 60), it % 2 ? 3 : 50);
    }

    check_rollback(50000, 50000, 1000, 300, true);
    return 0;
}
