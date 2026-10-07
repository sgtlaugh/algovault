#include "common.h"

#define main library_main
#include "../code_library/hashing.cpp"
#undef main

/// Hash equality against segment equality, collisions have probability ~n^2 / 2^61 per array and should never appear
int main(){
    /// modmul must return the canonical residue, a value in [mod, 2 mod) would hash equal segments differently
    for (long long it = 0; it < stress::scaled(1000000); it++){
        uint64_t a = stress::rng()() % mod, b = stress::rng()() % mod;
        if (it % 7 == 0) a = mod - 1 - stress::rand_int(0, 3), b = mod - 1 - stress::rand_int(0, 3);
        assert((uint64_t)modmul(a, b) == (uint64_t)((unsigned __int128)a * b % mod));
    }

    for (long long it = 0; it < stress::scaled(10000); it++){
        int n = stress::rand_int(1, it % 10 ? 30 : 2000), alphabet = stress::rand_int(1, 3);
        int mode = stress::rand_int(0, 2);
        vector<long long> a(n);
        for (auto& x : a){
            long long v = stress::rand_int(0, alphabet - 1);
            if (mode == 1) x = v - 1000000;                                    /// far below the old +997 offset
            else if (mode == 2) x = v ? (v == 1 ? LLONG_MIN : LLONG_MAX) : 0;  /// limits of the element type
            else x = v;
        }
        PolyHash h(a), hr(vector<long long>(a.rbegin(), a.rend()));

        for (int q = 0; q < 200; q++){
            int len = stress::rand_int(1, n), l1 = stress::rand_int(0, n - len), l2 = stress::rand_int(0, n - len);
            bool equal = std::equal(a.begin() + l1, a.begin() + l1 + len, a.begin() + l2);
            assert((h.get_hash(l1, l1 + len - 1) == h.get_hash(l2, l2 + len - 1)) == equal);
            assert(h.rev_hash(l1, l1 + len - 1) == hr.get_hash(n - l1 - len, n - l1 - 1));  /// the reversed segment in the reversed array
        }
    }
    return 0;
}
