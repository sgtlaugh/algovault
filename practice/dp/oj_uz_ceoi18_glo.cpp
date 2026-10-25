/***
 *
 * Longest Increasing Subsequence
 * LIS and LDS lengths ending at every index, plus one longest subsequence as indices
 *
 * Complexity: O(n log n)
 *
 * cmp is a strict weak ordering that defines "increasing", allow_equal = true gives the non-decreasing version
 * lis_vector(a)[i] = length of the longest increasing subsequence ending at index i
 * lis_indices(a) = increasing indices of one longest increasing subsequence, empty for empty input
 * lds_* are the decreasing counterparts, lis_indices(a, false, greater<T>()) gives one LDS
 *
 * Example:
 *   lis_vector(vector<int>{2, 8, 3, 9, 4})   // {1, 2, 2, 3, 3}
 *   lis_indices(vector<int>{2, 8, 3, 9, 4})  // {0, 2, 3}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
vector<int> lis_vector(const vector<T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    vector<T> idx;
    vector<int> res;

    for (const T &x : ar){
        auto it = allow_equal ? upper_bound(idx.begin(), idx.end(), x, cmp)
                              : lower_bound(idx.begin(), idx.end(), x, cmp);
        res.push_back(it - idx.begin() + 1);

        if (it == idx.end()) idx.push_back(x);
        else *it = x;
    }

    return res;
}

/// The rightmost i before j with len[i] == len[j] - 1 always precedes ar[j]: it was the last value written to
/// tail slot len[j] - 2 before j, and j landed in the slot after it, so no comparator check is needed
template <typename T, typename Compare = less<T>>
vector<int> lis_indices(const vector<T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    if (ar.empty()) return {};
    vector<int> len = lis_vector(ar, allow_equal, cmp);
    int last = max_element(len.begin(), len.end()) - len.begin();
    vector<int> res(len[last]);

    res.back() = last;
    for (int i = last - 1, k = len[last] - 2; k >= 0; i--){
        if (len[i] == k + 1) res[k--] = i;
    }

    return res;
}

template <typename T>
vector<int> lds_vector(const vector<T>& ar, bool allow_equal=false){
    return lis_vector(ar, allow_equal, [](const T& a, const T& b){ return b < a; });
}

template <typename T>
int lis_length(const vector<T>& ar, bool allow_equal=false){
    auto lis = lis_vector(ar, allow_equal);
    return lis.empty() ? 0 : *max_element(lis.begin(), lis.end());
}

template <typename T>
int lds_length(const vector<T>& ar, bool allow_equal=false){
    auto lds = lds_vector(ar, allow_equal);
    return lds.empty() ? 0 : *max_element(lds.begin(), lds.end());
}

int main(){
    int n;
    long long x;

    scanf("%d %lld", &n, &x);
    vector<long long> t(n);
    for (int i = 0; i < n; i++) scanf("%lld", &t[i]);

    vector<int> ending = lis_vector(t);
    vector<int> starting = lis_vector(vector<long long>(t.rbegin(), t.rend()), false, greater<long long>());
    reverse(starting.begin(), starting.end());

    vector<long long> vals = t;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());

    vector<int> fenwick(vals.size() + 1, 0);
    int res = 0;
    for (int j = 0; j < n; j++){
        int best = 0;
        for (int k = lower_bound(vals.begin(), vals.end(), t[j] + x) - vals.begin(); k > 0; k -= k & -k){
            best = max(best, fenwick[k]);
        }
        res = max(res, best + starting[j]);

        int pos = lower_bound(vals.begin(), vals.end(), t[j]) - vals.begin() + 1;
        for (; pos <= (int)vals.size(); pos += pos & -pos) fenwick[pos] = max(fenwick[pos], ending[j]);
    }

    printf("%d\n", res);
    return 0;
}
