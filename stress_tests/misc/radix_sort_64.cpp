#include "../common.h"

#define main library_main
#include "../../code_library/misc/radix_sort_64.cpp"
#undef main

int main(){
    vector<unsigned long long> buf;  /// reused across sizes, so a stale larger scratch is exercised
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 50 : 5000), byte = stress::rand_int(0, 7), mode = stress::rand_int(0, 3);
        vector<unsigned long long> v(n);
        for (auto& x : v){
            unsigned long long r = stress::rng()();
            /// Only one byte differs, every other pass must keep the order
            if (mode == 1) x = r & ((unsigned long long)255 << (8 * byte));
            else if (mode == 2) x = r % 4 ? 0 : ~(unsigned long long)0;      /// extremes with heavy duplication
            else if (mode == 3) x = r % 300;
            else x = r;
        }

        /// The same keys as signed long longs, shifted so that small values and the extremes straddle zero
        vector<long long> s(n);
        for (int i = 0; i < n; i++) s[i] = mode == 2 ? (v[i] ? LLONG_MAX : LLONG_MIN) : (long long)(v[i] - (mode == 3 ? 150 : 0));
        auto signed_expected = s;
        sort(signed_expected.begin(), signed_expected.end());
        radix_sort(s);
        assert(s == signed_expected);

        auto expected = v;
        sort(expected.begin(), expected.end());
        if (it % 3) radix_sort(v.data(), n, buf);  /// not it % 2, the large n iterations are all even
        else radix_sort(v.data(), n);
        assert(v == expected);
    }

    return 0;
}
