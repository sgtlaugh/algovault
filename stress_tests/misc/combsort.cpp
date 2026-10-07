#include "../common.h"

#define main library_main
#include "../../code_library/misc/combsort.c"
#undef main

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
        combsort(len, v.data());
        assert(v == expected);
    }
    return 0;
}
