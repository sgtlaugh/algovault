#include "../common.h"

#define main library_main
#include "../../code_library/dp/bounded_subset_sum.cpp"
#undef main

const int MAXW = 200000;

using Bits = bitset<MAXW + 1>;

/// Bitset DP adding one copy at a time, a value's copies stop once a copy reaches nothing new
Bits brute(const vector<pair<int, long long>>& items, int W){
    Bits can, mask;
    can[0] = 1;
    for (int s = 0; s <= W; s++) mask[s] = 1;

    for (const auto& [v, c] : items){
        for (long long k = 0; k < c; k++){
            Bits next = (can | (can << v)) & mask;
            if (next == can) break;
            can = next;
        }
    }

    return can;
}

void check_subset(const SubsetSum& ss, const vector<int>& a, int s){
    vector<int> ids = ss.subset(s);
    vector<char> used(a.size(), 0);
    long long sum = 0;
    for (int i : ids){
        assert(0 <= i && i < (int)a.size() && !used[i]);
        used[i] = 1, sum += a[i];
    }
    assert(sum == s);
}

void check(const vector<int>& a, int W, int subset_checks){
    vector<pair<int, long long>> items;
    for (int x : a) items.push_back({x, 1});
    Bits expected = brute(items, W);

    vector<char> can = bounded_subset_sums(items, W);
    assert((int)can.size() == W + 1);
    for (int s = 0; s <= W; s++) assert(can[s] == expected[s]);

    SubsetSum ss(a, W);
    for (int s = -2; s <= W + 2; s++) assert(ss.reachable(s) == (s >= 0 && s <= W && expected[s]));
    for (int s = 0; s <= W; s++){
        if (expected[s] && (subset_checks < 0 || stress::rand_int(0, W) < subset_checks)) check_subset(ss, a, s);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 14), W = stress::rand_int(0, it % 3 ? 60 : 400), hi = stress::rand_int(1, it % 2 ? 6 : W + 10);
        vector<int> a(n);
        for (auto& x : a) x = stress::rand_int(0, hi);
        check(a, W, -1);

        vector<pair<int, long long>> items(stress::rand_int(0, 6));
        for (auto& [v, c] : items){
            v = stress::rand_int(0, hi);
            c = stress::rand_int(0, 3) ? stress::rand_int(0, 6) : (long long)1e18;
        }
        vector<char> can = bounded_subset_sums(items, W);
        Bits expected = brute(items, W);
        for (int s = 0; s <= W; s++) assert(can[s] == expected[s]);
    }

    vector<int> a;
    for (int v = 1; v * (v + 1) / 2 <= MAXW; v++) a.push_back(v);
    shuffle(a.begin(), a.end(), stress::rng());
    check(a, MAXW, 300);

    a.clear();
    for (int left = MAXW; left > 0; ){
        int part = min<int>(left, stress::rand_int(1, 900));
        a.push_back(part), left -= part;
    }
    check(a, MAXW, 300);

    a.assign(5000, 0);
    for (auto& x : a) x = stress::rand_int(1, 40) * 1000 + stress::rand_int(0, 2);
    check(a, MAXW, 300);

    return 0;
}
