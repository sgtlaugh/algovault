/***
 *
 * Parallel Binary Search
 * For many queries at once, the fewest updates of a fixed sequence after which each query's condition holds
 *
 * Complexity: O(log(t + 2) (t U + q C + R)) for update cost U, check cost C and reset cost R, O(q + t) memory
 *
 * parallel_binary_search(q, t, reset, apply, check): res[i] = the fewest k in [0, t] whose updates make check(i) true
 *
 * Notes:
 *   - reset() restores the state before any update, apply(k) performs update k, check(i) reads query i's condition
 *   - Each round calls reset, then apply(0), apply(1), ... in order, checking every unfinished query once
 *   - res[i] = t + 1 when check(i) stays false after all t updates
 *   - check(i) must be monotone in the number of updates applied: once true it stays true
 *   - Use it when one binary search per query would replay the updates q log t times, as in SPOJ METEORS
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Reset, typename Apply, typename Check>
vector<int> parallel_binary_search(int q, int t, Reset reset, Apply apply, Check check){
    assert(q >= 0 && t >= 0);
    vector<int> lo(q, 0), hi(q, t + 1), head(t + 1), next(q);

    while (true){
        /// Queries are bucketed by their midpoint, so one pass over the updates answers all of them
        fill(head.begin(), head.end(), -1);
        int last = -1;
        for (int i = 0; i < q; i++){
            if (lo[i] == hi[i]) continue;
            int mid = (lo[i] + hi[i]) / 2;
            next[i] = head[mid], head[mid] = i, last = max(last, mid);
        }
        if (last == -1) break;

        reset();
        for (int k = 0; k <= last; k++){
            if (k > 0) apply(k - 1);
            for (int i = head[k]; i != -1; i = next[i]){
                if (check(i)) hi[i] = k;
                else lo[i] = k + 1;
            }
        }
    }

    return lo;
}

int main(){
    /// Updates add 3, 1, 4, 1, 5 to a running total; a query asks when the total first reaches its target
    vector<int> add = {3, 1, 4, 1, 5}, target = {0, 4, 9, 14, 15};
    long long total = 0;

    vector<int> res = parallel_binary_search(target.size(), add.size(),
        [&](){ total = 0; },
        [&](int k){ total += add[k]; },
        [&](int i){ return total >= target[i]; });

    assert(res[0] == 0);  /// 0 >= 0 before any update
    assert(res[1] == 2);  /// 3 + 1 = 4
    assert(res[2] == 4);  /// 3 + 1 + 4 + 1 = 9
    assert(res[3] == 5);  /// all five reach 14
    assert(res[4] == 6);  /// never, so t + 1
    return 0;
}
