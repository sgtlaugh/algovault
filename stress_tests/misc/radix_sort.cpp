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
            if (mode == 1) x = r & ((unsigned int)255 << (8 * byte));  /// only one byte differs, every other pass must keep the order
            else if (mode == 2) x = r % 4 ? 0 : ~(unsigned int)0;      /// extremes with heavy duplication
            else if (mode == 3) x = r % 300;
            else x = r;
        }

        auto expected = v;
        sort(expected.begin(), expected.end());
        if (it % 3) radix_sort(v.data(), n, buf);  /// not it % 2, the large n iterations are all even
        else radix_sort(v.data(), n);
        assert(v == expected);
    }

    return 0;
}
