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
    int n, q;
    for (int t = 1; scanf("%d %d", &n, &q) == 2 && (n || q); t++){
        vector<int> ar(n);
        for (int& x : ar) scanf("%d", &x);

        /// The longest decreasing run ending at i in the reversed array is the longest increasing one starting at i
        vector<int> start(ar.rbegin(), ar.rend());
        start = lis_vector(start, false, greater<int>());
        reverse(start.begin(), start.end());
        int len = *max_element(start.begin(), start.end());

        printf("Set %d:\n", t);
        for (int i = 1; i <= q; i++){
            int m;
            scanf("%d", &m);
            printf("  Subset %d:\n", i);
            if (m > len){
                puts("    Impossible");
                continue;
            }

            /// Greedy from the left: take the first index that can still start an increasing run of the needed length
            printf("   ");
            long long cur = LLONG_MIN;
            for (int j = 0; m; j++){
                if (start[j] >= m && ar[j] > cur) m--, cur = ar[j], printf(" %d", ar[j]);
            }
            puts("");
        }
    }
    return 0;
}
