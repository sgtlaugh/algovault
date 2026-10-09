/***
 *
 * Compress values in input vector in-place
 * Relative order of values are preserved if make_sorted=true
 *
 * 0 based indexing for compressed values
 *
 * Complexity: O(N) expected if sorting is not necessary (hash map), otherwise O(N log N) (no hashing)
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <class T>
void compress(vector<T>& v, bool make_sorted=true){
    if (!make_sorted){
        unordered_map <T, int> mp;
        mp.reserve(v.size());
        for (auto &&x: v) x = mp.emplace(x, mp.size()).first->second;
        return;
    }

    int n = v.size(), rank = 0;
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());
    for (int i = 0; i < n; i++){
        if (i && order[i - 1].first < order[i].first) rank++;
        v[order[i].second] = rank;
    }
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
