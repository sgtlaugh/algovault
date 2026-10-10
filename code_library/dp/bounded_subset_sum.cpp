/***
 *
 * Bounded Subset Sum
 * Which sums up to W are reachable by a sub-multiset of the items, with item ids for a chosen sum
 *
 * Complexity: bounded_subset_sums O(W * K / 64 + W), SubsetSum O(W * (D + 1) + n log n), both O(n + W) memory
 * D is the number of distinct values in [1, W], S the sum of all values with their counts
 * K <= 2 * min(W, D * (1 + log2 W), sqrt(2 * S)) is the number of copies bounded_subset_sums keeps after carrying
 * D distinct positive values sum to at least D^2 / 2, so D <= sqrt(2 * S), which makes both O(W sqrt W),
 * bounded_subset_sums with the / 64, when the values sum to O(W) (component sizes summing to n, a partition of n)
 *
 * bounded_subset_sums({{value, count}, ...}, W): can[s] for s in [0, W], counts up to 1e18, equal values merge
 * SubsetSum ss(a, W): ss.reachable(s) for any int s, ss.subset(s) returns distinct ids i with sum a[i] == s
 *     subset(s) requires reachable(s), subset(0) is empty
 *
 * Values and counts must be non-negative, W >= 0; values 0 and values > W are never used
 *
 * bounded_subset_sums keeps 1 or 2 copies of each v and carries every further pair to 2v, which reaches the same sums,
 * then ORs the bitset with itself shifted by v once per kept copy
 *
 * SubsetSum makes one pass per distinct value v, residue class r mod v at a time: walking s = r, r + v, r + 2v, ...
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
    vector<long long> cnt(W + 1, 0);
    for (const auto& [v, c] : items){
        assert(v >= 0 && c >= 0);
        if (v > 0 && v <= W) cnt[v] = min<long long>(W, cnt[v] + min<long long>(c, W));
    }

    int words = W / 64 + 1;
    vector<unsigned long long> bits(words, 0);
    bits[0] = 1;
    for (int v = 1; v <= W; v++){
        if (cnt[v] > 2){
            long long pairs = (cnt[v] - 1) / 2;
            cnt[v] -= 2 * pairs;
            if (v <= W / 2) cnt[2 * v] = min<long long>(W, cnt[2 * v] + pairs);
        }

        int q = v / 64, r = v % 64;
        for (int k = 0; k < cnt[v]; k++){
            if (r == 0){
                for (int i = words - 1; i >= q; i--) bits[i] |= bits[i - q];
            } else{
                for (int i = words - 1; i > q; i--) bits[i] |= bits[i - q] << r | bits[i - q - 1] >> (64 - r);
                bits[q] |= bits[0] << r;
            }
        }
    }

    vector<char> can(W + 1);
    for (int s = 0; s <= W; s++) can[s] = bits[s / 64] >> (s % 64) & 1;
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
    vector<char> can = bounded_subset_sums({{3, 2}, {5, 1}}, 12);  /// two 3s and one 5
    assert(can[0]);
    assert(can[11]);   /// 3 + 3 + 5
    assert(!can[9]);   /// would need three 3s

    SubsetSum ss({4, 6, 4, 9}, 15);
    assert(ss.reachable(14));
    assert(!ss.reachable(7));
    vector<int> ids = ss.subset(14);
    sort(ids.begin(), ids.end());
    assert((ids == vector<int>{0, 1, 2}));  /// 4 + 6 + 4, the ids come in any order
    return 0;
}
