#include "common.h"

#define main library_main
#include "../code_library/radix_sort_64.cpp"
#undef main

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 50 : 5000), byte = stress::rand_int(0, 7);
        vector<unsigned long long> v(n);
        for (auto& x : v){
            unsigned long long r = stress::rng()();
            int mode = it % 4;
            if (mode == 1) x = r & ((unsigned long long)255 << (8 * byte));  /// only one byte differs, every other pass must keep the order
            else if (mode == 2) x = r % 4 ? 0 : ~(unsigned long long)0;      /// extremes with heavy duplication
            else if (mode == 3) x = r % 300;
            else x = r;
        }

        auto expected = v;
        sort(expected.begin(), expected.end());
        radix_sort(v.data(), n);
        assert(v == expected);
    }
    return 0;
}
