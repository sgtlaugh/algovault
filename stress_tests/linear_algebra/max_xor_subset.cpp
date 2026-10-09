#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/max_xor_subset.cpp"
#undef main

/// Classic insertion basis indexed by leading bit, then a greedy descent
long long basis_max(const vector<long long>& ar){
    long long basis[64] = {0}, res = 0;
    for (long long x : ar){
        for (int b = 62; b >= 0 && x; b--){
            if (!(x >> b & 1)) continue;
            if (!basis[b]){
                basis[b] = x;
                x = 0;
            }
            else x ^= basis[b];
        }
    }
    for (int b = 62; b >= 0; b--) res = max(res, res ^ basis[b]);
    return res;
}

long long random_value(int bits){
    return bits ? (long long)(stress::rng()() >> (64 - bits)) : 0;
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(0, 14), bits = stress::rand_int(0, it % 3 ? 6 : 63);
        vector<long long> ar(n);
        for (auto& x : ar) x = random_value(stress::rand_int(0, bits));
        if (n && stress::rand_int(0, 3) == 0) ar[stress::rand_int(0, n - 1)] = LLONG_MAX;

        long long brute = 0;
        for (int mask = 0; mask < (1 << n); mask++){
            long long x = 0;
            for (int i = 0; i < n; i++) if (mask >> i & 1) x ^= ar[i];
            brute = max(brute, x);
        }
        assert(max_xor_subset(ar) == brute);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(15, 3000), bits = stress::rand_int(1, 63);
        vector<long long> ar(n);
        for (auto& x : ar) x = random_value(stress::rand_int(1, bits));
        assert(max_xor_subset(ar) == basis_max(ar));
    }
    return 0;
}
