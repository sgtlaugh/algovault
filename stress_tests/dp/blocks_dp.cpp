#include "../common.h"

#define main library_main
#include "../../code_library/dp/blocks_dp.c"
#undef main

/// Every click order: clicking removes one maximal run of equal colors
int brute(const vector<int>& v, map<vector<int>, int>& memo){
    if (v.empty()) return 0;
    auto found = memo.find(v);
    if (found != memo.end()) return found->second;

    int best = 0;
    for (size_t i = 0, j; i < v.size(); i = j){
        for (j = i; j < v.size() && v[j] == v[i]; j++){}
        vector<int> rest(v.begin(), v.begin() + i);
        rest.insert(rest.end(), v.begin() + j, v.end());
        best = max(best, (int)((j - i) * (j - i)) + brute(rest, memo));
    }

    return memo[v] = best;
}

int main(){
    /// The first block must skip its nearest equal neighbour to reach 50, random inputs this long rarely need that
    const int regression[] = {0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0};
    assert(get_max_score(12, regression) == 50);

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(0, it % 4 ? 9 : 12), colors = stress::rand_int(1, 3);  /// skipping a nearer equal block first pays off at n = 12
        vector<int> v(n);
        for (auto& x : v) x = stress::rand_int(1, colors) * (it % 2 ? 1 : 1000000);  /// colors are compared, never indexed
        map<vector<int>, int> memo;
        assert(get_max_score(n, v.data()) == brute(v, memo));
    }

    /// One color clears in a single click, distinct colors only one at a time, up to the 200 block limit
    for (long long it = 0; it <= stress::scaled(10); it++){
        int n = it ? stress::rand_int(1, 60) : 200;
        vector<int> same(n, 7), distinct(n);
        iota(distinct.begin(), distinct.end(), 0);
        assert(get_max_score(n, same.data()) == n * n);
        assert(get_max_score(n, distinct.data()) == n);
    }

    return 0;
}
