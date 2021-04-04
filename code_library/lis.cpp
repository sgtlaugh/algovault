/***
 * Fast O(n log n) algorithm to calculate LIS and LDS
 * LIS - Longest Increasing Subsequence
 * LDS - Longest Decreasing Subsequence
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// returns a vector L of size n, where n is the length of input vector
/// L[i] = length of longest increasing subsequence including and ending at index i
/// If allow_equal is true, longest non-decreasing subsequence is calculated instead

template <typename T>
vector <int> lis_vector(const vector <T>& ar, bool allow_equal=false){
    vector <T> idx;
    vector <int> res;

    for (const T &x : ar) {
        auto it = allow_equal ? upper_bound(idx.begin(), idx.end(), x)
                              : lower_bound(idx.begin(), idx.end(), x);
        res.push_back(it - idx.begin() + 1);

        if (it == idx.end()) idx.push_back(x);
        else *it = x;
    }

    return res;
}

/// returns a vector L of size n, where n is the length of input vector
/// L[i] = length of longest decreasing subsequence including and ending at index i
/// If allow_equal is true, longest non-increasing subsequence is calculated instead

template <typename T>
vector <int> lds_vector(const vector <T>& ar, bool allow_equal=true){
    int n = ar.size();

    vector <T> pos;
    for (int i = 0; i < n; i++) pos.push_back(ar[i]);
    sort(pos.begin(), pos.end());

    map <T, int> mp;
    for (int i = 0; i < n; i++) mp[pos[i]] = n - i - 1;

    vector <int> inv(n);
    for (int i = 0; i < n; i++) inv[i] = mp[ar[i]];
    return lis_vector(inv, allow_equal);
}

template <typename T>
int lis_length(const vector <T>& ar, bool allow_equal=false){
    auto lis = lis_vector(ar, allow_equal);
    return *max_element(lis.begin(), lis.end());
}

template <typename T>
int lds_length(const vector <T>& ar, bool allow_equal=false){
    auto lds = lds_vector(ar, allow_equal);
    return *max_element(lds.begin(), lds.end());
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

    return 0;
}
