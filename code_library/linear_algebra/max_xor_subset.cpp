/***
 *
 * Given an array of non-negative integers,
 * Finds a subset where the bitwise XOR of all the elements in the subset is maximum
 * Return the maximum xor value
 *
***/

#include <bits/stdc++.h>

using namespace std;

long long max_xor_subset(const vector<long long>& ar){
    vector<long long> basis;
    for (long long x : ar){
        for (long long b : basis) x = min(x, x ^ b);
        if (x) basis.push_back(x);
    }

    /// Each b lacks the leading bits of earlier ones, so no sort is needed, every leading bit ends up set
    long long res = 0;
    for (long long b : basis) res = max(res, res ^ b);
    return res;
}

int main(){
    assert(max_xor_subset({1, 2, 3}) == 3);
    assert(max_xor_subset({1, 2, 4}) == 7);
    assert(max_xor_subset({10, 4, 12, 23, 6, 60}) == 62);
    assert(max_xor_subset({10, 4, 12, 23, 97, 6, 6, 3, 51, 60}) == 127);

    return 0;
}
