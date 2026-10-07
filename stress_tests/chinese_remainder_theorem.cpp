#include "common.h"

#define main library_main
#include "../code_library/chinese_remainder_theorem.cpp"
#undef main

/// Pairwise coprime moduli whose product stays below 2^63, each below the documented 3 * 10^9
vector<int64_t> random_moduli(int count, int64_t max_mod){
    vector<int64_t> mods;
    __int128 prod = 1;
    for (int tries = 0; (int)mods.size() < count && tries < 100; tries++){
        int64_t m = stress::rand_int(1, max_mod);
        bool coprime = true;
        for (auto x : mods) coprime &= __gcd(x, m) == 1;
        if (coprime && prod * m < ((__int128)1 << 63)) mods.push_back(m), prod *= m;
    }
    return mods;
}

int main(){
    /// Tiny moduli: the unique answer in [0, prod) is found by brute force
    for (long long it = 0; it < stress::scaled(20000); it++){
        auto mods = random_moduli(stress::rand_int(1, 4), 12);
        int64_t prod = 1;
        for (auto m : mods) prod *= m;
        vector<int64_t> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-3 * m, 3 * m));

        auto solves = [&](int64_t x){
            for (size_t i = 0; i < mods.size(); i++) if (((x - rems[i]) % mods[i] + mods[i]) % mods[i]) return false;
            return true;
        };
        int64_t expected = 0;
        while (!solves(expected)) expected++;
        assert(CRT(rems, mods) == expected);
    }

    /// Large moduli up to the documented limits: check the congruences exactly
    for (long long it = 0; it < stress::scaled(100000); it++){
        auto mods = random_moduli(stress::rand_int(1, 3), it % 2 ? 3000000000LL : 1000000);
        __int128 prod = 1;
        for (auto m : mods) prod *= m;
        vector<int64_t> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-m, m));

        int64_t x = CRT(rems, mods);
        assert(0 <= x && x < prod);
        for (size_t i = 0; i < mods.size(); i++) assert((((__int128)x - rems[i]) % mods[i] + mods[i]) % mods[i] == 0);
    }
    return 0;
}
