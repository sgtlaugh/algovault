/***
 *
 * Given an array of positive integers,
 * Finds a subset where the bitwise XOR of all the elements in the subset is maximum
 * Return the maximum xor value
 *
***/

#include <bits/stdc++.h>

#define bitlen(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

using namespace std;

long long max_xor_subset(const vector<long long>& ar){
    long long m, x, res = 0;
    int i, j, l, n = ar.size();

    vector <long long> v[64];
    for (i = 0; i < n; i++) v[bitlen(ar[i])].push_back(ar[i]);

    for (i = 63; i > 0; i--){
        l = v[i].size();
        if (l){
            m = v[i][0];
            res = max(res, res ^ m);

            for (j = 1; j < l; j++){
                x = m ^ v[i][j];
                if (x) v[bitlen(x)].push_back(x);
            }
            v[i].clear();
        }
    }

    return res;
}

int main(){
    assert(max_xor_subset({1, 2, 3}) == 3);
    assert(max_xor_subset({1, 2, 4}) == 7);
    assert(max_xor_subset({10, 4, 12, 23, 6, 60}) == 62);
    assert(max_xor_subset({10, 4, 12, 23, 97, 6, 6, 3, 51, 60}) == 127);
    return 0;
}
