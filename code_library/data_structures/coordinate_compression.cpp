/***
 *
 * Coordinate Compression
 * Compress values in input vector in-place
 * Relative order of values are preserved if make_sorted=true, otherwise ids follow first appearance
 *
 * 0 based indexing for compressed values
 * Returns the distinct values in id order, vals[v[i]] is the original value of v[i]
 *
 * Complexity: O(N) expected if sorting is not necessary (hash map), otherwise O(N log N) (no hashing)
 *
 * Example:
 *   auto vals = compress(v);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
vector<T> compress(vector<T>& v, bool make_sorted=true){
    vector<T> vals;
    if (!make_sorted){
        unordered_map<T, int> mp;
        mp.reserve(v.size());
        for (auto &&x: v){
            auto [it, fresh] = mp.emplace(x, vals.size());
            if (fresh) vals.push_back(x);
            x = it->second;
        }
        return vals;
    }

    int n = v.size();
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());
    for (int i = 0; i < n; i++){
        if (!i || order[i - 1].first < order[i].first) vals.push_back(order[i].first);
        v[order[i].second] = vals.size() - 1;
    }
    return vals;
}

int main(){
    vector<int> u = {2000000000, 1000000000, 2000000000, 1, 10, 5};
    auto u_vals = compress(u);
    assert(u == vector<int>({4, 3, 4, 0, 2, 1}));
    assert(u_vals == vector<int>({1, 5, 10, 1000000000, 2000000000}));

    vector<long long> v = {2000000000000LL, 1000000000, 2000000000000LL, 1, 10, 5};
    auto v_vals = compress(v, false);
    assert(v == vector<long long>({0, 1, 0, 2, 3, 4}));
    assert(v_vals == vector<long long>({2000000000000LL, 1000000000, 1, 10, 5}));

    vector<int> empty;
    assert(compress(empty).empty() && compress(empty, false).empty());

    return 0;
}
