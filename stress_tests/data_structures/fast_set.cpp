#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fast_set.cpp"
#undef main

void check_queries(const FastSet& fs, const set<int>& ref, int x){
    int n = fs.n;
    auto it = ref.lower_bound(x);
    assert(fs.next(x) == (it == ref.end() ? n : *it));

    it = ref.upper_bound(x);
    assert(fs.prev(x) == (it == ref.begin() ? -1 : *prev(it)));
    assert(fs.contains(x) == (ref.count(x) == 1));
}

/// Random inserts and erases, dense or sparse, checked at random points and at word and level edges
void random_ops(int n, int ops, int density, FastSet::Fill fill){
    FastSet fs(n, fill);
    set<int> ref;
    for (int x = 0; fill == FastSet::FULL && x < n; x++) ref.insert(x);
    vector<int> edges = {INT_MIN, -65, -1, 0, n - 1, n, n + 1, INT_MAX};
    for (int b = 64; b <= n; b *= 64){
        for (int k = b; k < n && k <= 4 * b; k += b) edges.insert(edges.end(), {k - 1, k, k + 1});
    }

    for (int it = 0; it < ops; it++){
        int x = n ? stress::rand_int(0, n - 1) : 0;
        if (n && stress::rand_int(0, 99) < density) fs.insert(x), ref.insert(x);
        else if (n) fs.erase(x), ref.erase(x);

        check_queries(fs, ref, stress::rand_int(-2, n + 1));
        check_queries(fs, ref, edges[stress::rand_int(0, edges.size() - 1)]);
        if (!ref.empty()){
            int member = *ref.lower_bound(stress::rand_int(0, *ref.rbegin()));
            check_queries(fs, ref, member - 1);
            check_queries(fs, ref, member + 1);
        }
    }
}

int main(){
    for (FastSet empty : {FastSet(0), FastSet(0, FastSet::FULL)}){
        assert(empty.next(0) == 0);
        assert(empty.prev(5) == -1);
        assert(!empty.contains(0));
    }

    /// Every subset of [0, n) for n <= 10, every query point, built by inserts or by erasing from FULL
    for (int n = 1; n <= 10; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            FastSet built(n), carved(n, FastSet::FULL);
            set<int> ref;
            for (int x = 0; x < n; x++){
                if (mask >> x & 1) built.insert(x), ref.insert(x);
            }
            for (int x = 0; x < n; x++){
                if (!(mask >> x & 1)) carved.erase(x);
            }
            for (int x = -2; x <= n + 1; x++) check_queries(built, ref, x), check_queries(carved, ref, x);
        }
    }

    /// Sizes around word and level boundaries, 64 = one full word, 4096 = one full second level
    for (int n : {1, 2, 63, 64, 65, 127, 128, 129, 4095, 4096, 4097, 262143, 262144, 262145}){
        for (int density : {2, 50, 98}) random_ops(n, stress::scaled(3000), density, FastSet::EMPTY);
        random_ops(n, stress::scaled(3000), 10, FastSet::FULL);
    }
    for (int it = 0; it < stress::scaled(200); it++){
        FastSet::Fill fill = it % 2 ? FastSet::FULL : FastSet::EMPTY;
        random_ops(stress::rand_int(1, 300), 200, stress::rand_int(0, 100), fill);
    }

    /// A 4-level tree, n = 1e7: lone members far apart force climbs to the top and back down
    int n = 10000000;
    FastSet big(n);
    set<int> ref;
    for (int x : {0, 4095, 4096, 262144, 5000000, n - 1}) big.insert(x), ref.insert(x);
    for (int it = 0; it < stress::scaled(20000); it++){
        int x = stress::rand_int(0, n - 1);
        if (it % 3 == 0) big.insert(x), ref.insert(x);
        auto member = ref.lower_bound(x);
        if (it % 3 == 1 && member != ref.end()) big.erase(*member), ref.erase(member);
        check_queries(big, ref, stress::rand_int(-1, n));
    }

    /// FULL at n = 1e7 with the last 64 members erased: a stray bit past n in any level makes next leave the tree
    FastSet full(n, FastSet::FULL);
    for (int x = n - 64; x < n; x++) full.erase(x);
    assert(full.next(n - 100) == n - 100);
    assert(full.next(n - 64) == n);
    assert(full.prev(n + 5) == n - 65);

    stress::expect_assert("data_structures/fast_set.cpp", []{ FastSet bad(64); bad.insert(64); });
    stress::expect_assert("data_structures/fast_set.cpp", []{ FastSet bad(64); bad.erase(-1); });

    return 0;
}
