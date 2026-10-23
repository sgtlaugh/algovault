/***
 *
 * Coordinate Compression
 * Replaces each value of a vector in place with a small id, and returns the distinct values
 *
 * Complexity: O(n log n), no hashing
 *
 * compress(v): ids follow value order, returns the distinct values sorted, so vals[v[i]] is the original v[i]
 * compress_by_appearance(v): ids follow first appearance, returns the distinct values in that order
 *
 * Ids run from 0 to the number of distinct values minus 1
 * The element type must hold ids up to n - 1, 65536 distinct values overflow a vector<short>
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
vector<T> compress(vector<T>& v){
    int n = v.size();
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());

    vector<T> vals;
    vector<int> ids(n);
    for (int i = 0; i < n; i++){
        if (!i || order[i - 1].first < order[i].first) vals.push_back(order[i].first);
        ids[order[i].second] = vals.size() - 1;
    }

    for (int i = 0; i < n; i++) v[i] = ids[i];
    return vals;
}

template <class T>
vector<T> compress_by_appearance(vector<T>& v){
    vector<T> sorted_vals = compress(v), vals;
    vector<int> relabel(sorted_vals.size(), -1);

    /// Relabel by scanning in input order instead of hashing, std::hash is the identity on ints and hackable
    for (T& x : v){
        int id = x;
        if (relabel[id] == -1){
            relabel[id] = vals.size();
            vals.push_back(sorted_vals[id]);
        }
        x = relabel[id];
    }
    return vals;
}

int main(){
    /// Ids follow value order
    vector<int> u = {2000000000, 1000000000, 2000000000, 1, 10, 5};
    vector<int> u_vals = compress(u);
    assert(u == vector<int>({4, 3, 4, 0, 2, 1}));
    assert(u_vals == vector<int>({1, 5, 10, 1000000000, 2000000000}));

    /// Ids follow first appearance
    vector<long long> v = {2000000000000LL, 1000000000, 2000000000000LL, 1, 10, 5};
    vector<long long> v_vals = compress_by_appearance(v);
    assert(v == vector<long long>({0, 1, 0, 2, 3, 4}));
    assert(v_vals == vector<long long>({2000000000000LL, 1000000000, 1, 10, 5}));

    return 0;
}
