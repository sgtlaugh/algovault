#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/persistent_segment_tree.cpp"
#undef main

/// Every version against a full copy of the array, branching updates from random old versions
template<typename Merge>
void check_versions(int n, int ops, long long range, bool sum){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-range, range);
    PersistentSegmentTree<long long, Merge> tree(a);
    vector<vector<long long>> copies = {a};

    for (int q = 0; q < ops; q++){
        int v = stress::rand_int(0, 2) ? copies.size() - 1 : stress::rand_int(0, copies.size() - 1);
        int type = stress::rand_int(0, 3), i = stress::rand_int(0, n - 1);
        if (type <= 1){
            vector<long long> next = copies[v];
            long long x = stress::rand_int(-range, range);
            if (type == 0 && sum) next[i] += x, assert(tree.add(v, i, x) == (int)copies.size());
            else next[i] = x, assert(tree.set(v, i, x) == (int)copies.size());
            copies.push_back(next);
            continue;
        }

        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        long long expected = copies[v][l];
        for (int j = l + 1; j <= r; j++) expected = Merge()(expected, copies[v][j]);
        assert(tree.query(v, l, r) == expected);
        assert(tree.get(v, i) == copies[v][i]);
    }

    for (int v = 0; v < (int)copies.size(); v += 1 + (int)copies.size() / 20){
        for (int i = 0; i < n; i++) assert(tree.get(v, i) == copies[v][i]);
    }
}

/// Persistent array use with values at the long long extremes, which a summing Merge would overflow
void check_array(int n, int ops){
    vector<long long> pool = {LLONG_MIN, LLONG_MAX, LLONG_MIN + 1, LLONG_MAX - 1, -1, 0, 1};
    vector<long long> a(n);
    for (auto& x : a) x = pool[stress::rand_int(0, pool.size() - 1)];
    PersistentSegmentTree<long long, FirstMerge> arr(a);
    vector<vector<long long>> copies = {a};

    for (int q = 0; q < ops; q++){
        int v = stress::rand_int(0, copies.size() - 1), i = stress::rand_int(0, n - 1);
        if (stress::rand_int(0, 1)){
            vector<long long> next = copies[v];
            next[i] = pool[stress::rand_int(0, pool.size() - 1)];
            arr.set(v, i, next[i]);
            copies.push_back(next);
        }
        else assert(arr.get(v, i) == copies[v][i]);
    }
}

/// Direct kth on any pair lo <= hi of a linear history of nonnegative adds, against a prefix scan of the difference
void check_kth_versions(int n, int ops){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(0, 3);
    PersistentSegmentTree<long long> tree(a);
    vector<vector<long long>> copies = {a};

    for (int q = 0; q < ops; q++){
        vector<long long> next = copies.back();
        int i = stress::rand_int(0, n - 1);
        long long x = stress::rand_int(0, 3);
        next[i] += x, tree.add(copies.size() - 1, i, x);
        copies.push_back(next);

        int lo = stress::rand_int(0, copies.size() - 1), hi = stress::rand_int(lo, copies.size() - 1);
        long long total = 0;
        for (int j = 0; j < n; j++) total += copies[hi][j] - copies[lo][j];
        if (total == 0) continue;

        long long k = stress::rand_int(0, total - 1), prefix = 0;
        int expected = 0;
        while ((prefix += copies[hi][expected] - copies[lo][expected]) <= k) expected++;
        assert(tree.kth(lo, hi, k) == expected);
    }
}

/// Order statistics against sorting the range directly
void check_kth(int n, int queries, long long range){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-range, range);
    if (n > 2 && stress::rand_int(0, 3) == 0) a[0] = LLONG_MIN, a[n - 1] = LLONG_MAX;
    RangeKth<long long> rk(a);

    for (int q = 0; q < queries; q++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (q % 10 == 0) l = 0, r = n - 1;
        vector<long long> part(a.begin() + l, a.begin() + r + 1);
        sort(part.begin(), part.end());
        int k = stress::rand_int(0, r - l);
        assert(rk.kth(l, r, k) == part[k]);
    }
}

int main(){
    /// Empty inputs still build version 0
    assert(PersistentSegmentTree<int>({}).versions() == 1);
    assert(RangeKth<int>({}).counts.versions() == 1);

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it % 40 + 1 : stress::rand_int(1, 200);
        long long range = it % 3 ? 5 : 1000000000000LL;
        check_versions<plus<long long>>(n, 200, range, true);
        check_versions<MinMerge>(n, 200, range, false);
        check_versions<FirstMerge>(n, 200, range, false);
        check_array(n, 200);
        check_kth_versions(n, 100);
        check_kth(n, 60, it % 2 ? 5 : 1000000000000000000LL);
    }

    for (int n : {1, 2, 3, 127, 128, 129}) check_versions<plus<long long>>(n, 3000, 3, true), check_kth(n, 2000, 3);
    check_versions<plus<long long>>(20000, 400, 1000000000, true);
    check_kth(100000, 300, 1000000000);

    return 0;
}
