#include <bits/stdtr1c++.h>

#define bitlen(x) ((x) == 0 ? (0) : (64 - __builtin_clzll(x)))

using namespace std;

/***
 *
 * given an array of positive integers, find a subset where the bitwise XOR of all the elements in the subset is maximum
 * return the maximum xor
 *
***/

long long max_xor_subset(int n, long long* ar){
    vector <long long> v[64];
    for (int i = 0; i < n; i++) v[bitlen(ar[i])].push_back(ar[i]);

    long long m, x, res = 0;
    for (int i = 63; i > 0; i--){
        int l = v[i].size();
        if (l){
            m = v[i][0];
            res = max(res, res ^ m);

            for (int j = 1; j < l; j++){
                x = m ^ v[i][j];
                if (x) v[bitlen(x)].push_back(x);
            }
            v[i].clear();
        }
    }

    return res;
}

int main(){
    long long A[] = {1, 2, 3};
    printf("%lld\n", max_xor_subset(3, A)); /// 3 (take only 3)

    long long B[] = {1, 2, 4};
    printf("%lld\n", max_xor_subset(3, B)); /// 7 (take all - 1, 2 and 4)

    return 0;
}
