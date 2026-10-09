/***
 *
 * Bounded Subset Sum
 * Which sums up to W are reachable by a sub-multiset of the items, with item ids for a chosen sum
 *
 * Complexity: O(W * (D + 1) + n log n) time, O(n + W) memory, where D is the number of distinct values in [1, W]
 * D distinct positive values sum to at least D^2 / 2, so D <= sqrt(2 * S) for S the sum of all values,
 * which makes it O(W sqrt W) when the values sum to O(W) (component sizes summing to n, a partition of n)
 *
 * bounded_subset_sums({{value, count}, ...}, W): can[s] for s in [0, W], counts up to 1e18, equal values merge
 * SubsetSum ss(a, W): ss.reachable(s) for any int s, ss.subset(s) returns distinct ids i with sum a[i] == s
 *     subset(s) requires reachable(s), subset(0) is empty
 *
 * Values and counts must be non-negative, W >= 0; values 0 and values > W are never used
 *
 * One pass per distinct value v, residue class r mod v at a time: walking s = r, r + v, r + 2v, ...
 * every already reachable s refills the budget to count(v), and each unreachable s spends one copy.
 * The newly reached s keep the id of the copy they spent, so subset(s) walks back to an older sum
 *
 * Example:
 *   SubsetSum ss({4, 6, 4, 9}, 15);
 *   ss.reachable(14), ss.subset(14)   // true, ids {0, 1, 2} in some order
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<char> bounded_subset_sums(const vector<pair<int, long long>>& items, int W){
    assert(W >= 0);
    vector<int> cnt(W + 1, 0);
    for (const auto& [v, c] : items){
        assert(v >= 0 && c >= 0);
        if (v > 0 && v <= W) cnt[v] = min<long long>(W, cnt[v] + min<long long>(c, W));
    }

    vector<char> can(W + 1, 0);
    can[0] = 1;
    for (int v = 1; v <= W; v++){
        if (!cnt[v]) continue;
        for (int r = 0; r < v; r++){
            int left = 0;
            for (int s = r; s <= W; s += v){
                if (can[s]) left = cnt[v];
                else if (left > 0) can[s] = 1, left--;
            }
        }
    }

    return can;
}

struct SubsetSum{
    vector<int> a, last;
    vector<char> can;

    SubsetSum(const vector<int>& a, int W) : a(a){
        assert(W >= 0);
        last.assign(W + 1, -1);
        can.assign(W + 1, 0);

        vector<int> order;
        for (int i = 0; i < (int)a.size(); i++){
            assert(a[i] >= 0);
            if (a[i] > 0 && a[i] <= W) order.push_back(i);
        }
        sort(order.begin(), order.end(), [&](int i, int j){ return a[i] < a[j]; });

        can[0] = 1;
        for (int lo = 0, hi; lo < (int)order.size(); lo = hi){
            int v = a[order[lo]];
            for (hi = lo; hi < (int)order.size() && a[order[hi]] == v; hi++) {}

            for (int r = 0; r < v; r++){
                int left = 0;
                for (int s = r; s <= W; s += v){
                    if (can[s]) left = hi - lo;
                    else if (left > 0) can[s] = 1, last[s] = order[lo + --left];
                }
            }
        }
    }

    bool reachable(int s) const{
        return s >= 0 && s < (int)can.size() && can[s];
    }

    vector<int> subset(int s) const{
        assert(reachable(s));
        vector<int> ids;
        for (; s > 0; s -= a[ids.back()]) ids.push_back(last[s]);
        return ids;
    }
};

int main(){
    auto as_string = [](const vector<char>& can){
        string res;
        for (char c : can) res += c ? '1' : '0';
        return res;
    };

    auto sorted_subset = [](const SubsetSum& ss, int s){
        vector<int> ids = ss.subset(s);
        sort(ids.begin(), ids.end());
        return ids;
    };

    assert(as_string(bounded_subset_sums({{2, 1}, {3, 1}, {6, 1}, {7, 1}}, 10)) == "10110111111");
    assert(as_string(bounded_subset_sums({{3, 2}, {5, 1}}, 12)) == "1001011010010");
    assert(as_string(bounded_subset_sums({{2, 1}, {2, 1}}, 5)) == "101010");
    assert(as_string(bounded_subset_sums({{4, (long long)1e18}}, 10)) == "10001000100");
    assert(as_string(bounded_subset_sums({{1, (long long)1e18}, {1, (long long)1e18}}, 5)) == "111111");
    assert(as_string(bounded_subset_sums({{0, 5}, {11, 3}, {7, 0}}, 3)) == "1000");
    assert(as_string(bounded_subset_sums({}, 0)) == "1");

    SubsetSum ss({4, 6, 4, 9}, 15);
    string reach;
    for (int s = 0; s <= 15; s++) reach += ss.reachable(s) ? '1' : '0';
    assert(reach == "1000101011100111");
    assert(!ss.reachable(-1) && !ss.reachable(16) && !ss.reachable(17));
    assert((sorted_subset(ss, 14) == vector<int>{0, 1, 2}));
    assert((sorted_subset(ss, 15) == vector<int>{1, 3}));
    assert((sorted_subset(ss, 8) == vector<int>{0, 2}));
    assert(ss.subset(0).empty());

    SubsetSum skip({0, 7, 20, 3}, 10);
    assert((sorted_subset(skip, 10) == vector<int>{1, 3}));
    assert(!skip.reachable(4) && !skip.reachable(20));

    SubsetSum none({}, 0);
    assert(none.reachable(0) && none.subset(0).empty() && !none.reachable(1));

    return 0;
}
