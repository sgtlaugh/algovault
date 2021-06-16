/***
 *
 * Compress values in input vector in-place
 * Relative order of values are preserved if make_sorted=true
 *
 * 0 based indexing for compressed values
 *
 * Complexity: O(N) if sorting is not necessary, otherwise O(N log N)
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <class T>
void compress(vector<T>& v, bool make_sorted=true){
    unordered_map <T, int> mp;

    if (!make_sorted){
        for (auto &&x: v) x = mp.emplace(x, mp.size()).first->second;
        return;
    }

    vector<T> u = v;
    sort(u.begin(), u.end());
    for (auto &&x: u) mp.emplace(x, mp.size()).first->second;
    for (auto &&x: v) x = mp[x];
}

int main(){
    vector <int> u = {2000000000, 1000000000, 2000000000, 1, 10, 5};
    compress(u);
    assert(u == vector<int>({4, 3, 4, 0, 2, 1}));

    vector <long long> v = {2000000000000LL, 1000000000, 2000000000000LL, 1, 10, 5};
    compress(v, false);
    assert(v == vector<long long>({0, 1, 0, 2, 3, 4}));

    return 0;
}
