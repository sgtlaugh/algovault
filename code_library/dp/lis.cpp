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
vector <int> lis_vector(const vector <T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    vector <T> idx;
    vector <int> res;

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
vector <int> lis_indices(const vector <T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    if (ar.empty()) return {};
    vector <int> len = lis_vector(ar, allow_equal, cmp);
    int last = max_element(len.begin(), len.end()) - len.begin();
    vector <int> res(len[last]);

    res.back() = last;
    for (int i = last - 1, k = len[last] - 2; k >= 0; i--){
        if (len[i] == k + 1) res[k--] = i;
    }

    return res;
}

template <typename T>
vector <int> lds_vector(const vector <T>& ar, bool allow_equal=false){
    return lis_vector(ar, allow_equal, [](const T& a, const T& b){ return b < a; });
}

template <typename T>
int lis_length(const vector <T>& ar, bool allow_equal=false){
    auto lis = lis_vector(ar, allow_equal);
    return lis.empty() ? 0 : *max_element(lis.begin(), lis.end());
}

template <typename T>
int lds_length(const vector <T>& ar, bool allow_equal=false){
    auto lds = lds_vector(ar, allow_equal);
    return lds.empty() ? 0 : *max_element(lds.begin(), lds.end());
}

int main(){
    vector <int> ar, res;

    ar = {1, 2, 4, 3};
    res = {1, 2, 3, 3};
    assert(lis_vector(ar) == res);
    assert(lis_length(ar) == 3);

    ar = {4, 3, 5, 2, 1};
    res = {1, 2, 1, 3, 4};
    assert(lds_vector(ar) == res);
    assert(lds_length(ar) == 4);

    ar = {4, 3, 5, 2, 1, 5, 5, 5, 5, 4};
    res = {1, 2, 1, 3, 4, 2, 3, 4, 5, 6};
    assert(lds_vector(ar, true) == res);
    assert(lds_length(ar, true) == 6);

    ar = {3, 1, 4, 4, 6, 7, 8, 10, 1, 2, 2, 2, 2, 3, 8, 9};
    res = {1, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 8, 9};
    assert(lis_vector(ar, true) == res);
    assert(lis_length(ar, true) == 9);

    assert((lis_vector(vector<int>{2, 8, 3, 9, 4}) == vector<int>{1, 2, 2, 3, 3}));
    assert((lis_indices(vector<int>{2, 8, 3, 9, 4}) == vector<int>{0, 2, 3}));
    assert((lis_indices(vector<int>{}) == vector<int>{}));
    assert((lis_indices(vector<int>{-7}) == vector<int>{0}));
    assert((lis_indices(vector<int>{5, 5, 5}) == vector<int>{0}));
    assert((lis_indices(vector<int>{5, 5, 5}, true) == vector<int>{0, 1, 2}));
    assert((lis_indices(vector<int>{2, 8, 3, 9, 4, 1, 5}) == vector<int>{0, 2, 4, 6}));
    assert((lis_indices(vector<int>{1, 3, 3, 2, 2, 2, 4}, true) == vector<int>{0, 3, 4, 5, 6}));
    assert((lis_indices(vector<int>{9, 4, 10, 3, 1}, false, greater<int>()) == vector<int>{0, 1, 3, 4}));
    assert((lis_indices(vector<string>{"ccc", "a", "bb", "dddd"}, false, [](const string& x, const string& y){ return x.size() < y.size(); }) == vector<int>{1, 2, 3}));

    return 0;
}
