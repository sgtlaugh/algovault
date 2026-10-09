#include "../common.h"

#define main library_main
#include "../../code_library/misc/combsort.cpp"
#undef main

/// combsort on vectors, raw arrays and long long keys against std::sort
int main(){
    for (long long it = 0; it < stress::scaled(8000); it++){
        int len = stress::rand_int(0, it % 10 ? 40 : 3000), mode = it % 4;
        vector<int> v(len);
        for (int i = 0; i < len; i++){
            if (mode == 0) v[i] = stress::rand_int(INT_MIN, INT_MAX);
            else if (mode == 1) v[i] = stress::rand_int(0, 3);
            else if (mode == 2) v[i] = len - i;                    /// reversed, the worst case for small gaps
            else v[i] = i + (stress::rand_int(0, 9) ? 0 : stress::rand_int(-5, 5));
        }

        auto expected = v;
        sort(expected.begin(), expected.end());
        auto raw = v;
        combsort(v.begin(), v.end());
        combsort(raw.data(), raw.data() + len);
        assert(v == expected);
        assert(raw == expected);

        vector<long long> wide(stress::rand_int(0, 60));
        for (auto& x : wide) x = stress::rand_int(LLONG_MIN, LLONG_MAX);
        auto wide_expected = wide;
        sort(wide_expected.begin(), wide_expected.end());
        combsort(wide.begin(), wide.end());
        assert(wide == wide_expected);
    }

    for (int n = 0; n <= 8; n++){
        vector<int> perm(n);
        iota(perm.begin(), perm.end(), 0);
        do {
            auto v = perm;
            combsort(v.begin(), v.end());
            for (int i = 0; i < n; i++) assert(v[i] == i);
        } while (next_permutation(perm.begin(), perm.end()));
    }

    /// 1e6 random values, a gap sequence that degrades to bubble sort would not finish
    vector<int> big(1000000);
    for (auto& x : big) x = stress::rand_int(0, 1000000006);
    auto expected = big;
    sort(expected.begin(), expected.end());
    combsort(big.begin(), big.end());
    assert(big == expected);

    return 0;
}
