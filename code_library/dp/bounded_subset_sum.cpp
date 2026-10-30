/***
 *
 * Bounded Subset Sum
 * Which sums up to W are reachable by a sub-multiset of the items, with item ids for a chosen sum
 *
 * Complexity:
 *   - bounded_subset_sums: O(W K / 64 + W) time, O(n + W) memory, K <= 2 min(W, D (1 + log2 W), sqrt(2 S)) copies
 *   - SubsetSum: O(W (D + 1) + n log n) time, O(n + W) memory
 *
 * bounded_subset_sums(items, W): can[s] for s in [0, W], items any range of (value, count) pairs like map<int, int>
 * SubsetSum ss(a, W): sub-multisets of the values a[i] up to W
 * ss.reachable(s): whether some sub-multiset sums to s, for any int s
 * ss.subset(s): distinct ids i with the a[i] summing to s, needs reachable(s), empty for s = 0
 *
 * Notes:
 *   - D is the number of distinct values in [1, W], S the sum of all values with their counts
 *   - D <= sqrt(2 S), so both are O(W sqrt W) when the values sum to O(W), like component sizes summing to n
 *   - Values and counts are non-negative, counts up to 1e18, W >= 0, and equal values add their counts
 *   - Values 0 and values above W are never used
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Items = vector<pair<int, long long>>>
vector<char> bounded_subset_sums(const Items& items, int W){
    assert(W >= 0);
    vector<long long> cnt(W + 1, 0);
    for (const auto& item : items){
        long long v = item.first, c = item.second;
        assert(v >= 0 && c >= 0);
        if (v > 0 && v <= W) cnt[v] = min<long long>(W, cnt[v] + min<long long>(c, W));
    }

    int words = W / 64 + 1, top = 0;  /// words above top are still all zero
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
            int hi = min(words - 1, top + q + 1);
            if (r == 0){
                for (int i = hi; i >= q; i--) bits[i] |= bits[i - q];
            } else{
                for (int i = hi; i > q; i--) bits[i] |= bits[i - q] << r | bits[i - q - 1] >> (64 - r);
                bits[q] |= bits[0] << r;
            }
            top = hi;
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

    map<int, int> freq = {{2, 3}, {7, 1}};
    assert((bounded_subset_sums(freq, 13) == vector<char>{1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1}));

    SubsetSum ss({4, 6, 4, 9}, 15);
    assert(ss.reachable(14));
    assert(!ss.reachable(7));
    vector<int> ids = ss.subset(14);
    sort(ids.begin(), ids.end());
    assert((ids == vector<int>{0, 1, 2}));  /// 4 + 6 + 4, the ids come in any order
    return 0;
}
