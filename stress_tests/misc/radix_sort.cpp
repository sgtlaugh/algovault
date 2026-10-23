#include "../common.h"

#define main library_main
#include "../../code_library/misc/radix_sort.cpp"
#undef main

int main(){
    vector<unsigned int> buf;  /// reused across sizes, so a stale larger scratch is exercised
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 50 : 5000), byte = stress::rand_int(0, 3), mode = stress::rand_int(0, 3);
        vector<unsigned int> v(n);
        for (auto& x : v){
            unsigned int r = stress::rng()();
            /// Only one byte differs, every other pass must keep the order
            if (mode == 1) x = r & ((unsigned int)255 << (8 * byte));
            else if (mode == 2) x = r % 4 ? 0 : ~(unsigned int)0;      /// extremes with heavy duplication
            else if (mode == 3) x = r % 300;
            else x = r;
        }

        /// The same keys as signed ints, shifted so that small values and the extremes straddle zero
        vector<int> s(n);
        for (int i = 0; i < n; i++) s[i] = mode == 2 ? (v[i] ? INT_MAX : INT_MIN) : (int)(v[i] - (mode == 3 ? 150 : 0));
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
