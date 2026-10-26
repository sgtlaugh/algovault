/***
 *
 * Longest Increasing Subsequence
 * LIS and LDS lengths ending at every index, plus one longest subsequence as indices
 *
 * Complexity: O(n log n)
 *
 * lis_vector(a, allow_equal, cmp)[i]: length of the longest increasing subsequence ending at index i
 * lis_indices(a, allow_equal, cmp): increasing indices of one longest increasing subsequence, empty for empty input
 * lis_length(a, allow_equal): length of the longest increasing subsequence
 * lds_vector(a, allow_equal)[i]: length of the longest decreasing subsequence ending at index i
 * lds_length(a, allow_equal): length of the longest decreasing subsequence
 *
 * Notes:
 *   - cmp is a strict weak ordering that defines "increasing", less<T> by default
 *   - allow_equal = true gives the non-decreasing and non-increasing versions
 *   - lis_indices(a, false, greater<T>()) gives one longest decreasing subsequence
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
vector<int> lis_vector(const vector<T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    vector<T> tails;
    vector<int> res;
    auto before = [&](const T& tail, const T& x){ return allow_equal ? !cmp(x, tail) : cmp(tail, x); };

    for (const T &x : ar){
        if (tails.empty() || before(tails.back(), x)){
            tails.push_back(x);
            res.push_back(tails.size());
            continue;
        }

        /// Branchless search in [0, size - 1], the back is not before x; 2.5x faster than lower_bound at random n = 1e6
        const T* first = tails.data();
        size_t n = tails.size();
        while (n > 1){
            size_t half = n / 2;
            first += half * before(first[half - 1], x);
            n -= half;
        }

        res.push_back(first - tails.data() + 1);
        tails[first - tails.data()] = x;
    }

    return res;
}

/// The rightmost i before j with len[i] == len[j] - 1 always precedes ar[j], so no comparator check is needed
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
    vector<int> ar, res;

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
