#include "../common.h"

#define main library_main
#include "../../code_library/algebra/walsh_hadamard.cpp"
#undef main

int main(){
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = 1 << stress::rand_int(0, it % 20 ? 7 : 11), range = stress::rand_int(0, 1) ? 3 : 100000;
        vector<long long> a(n), b(n);
        for (auto& x : a) x = stress::rand_int(-range, range);
        for (auto& x : b) x = stress::rand_int(-range, range);

        vector<long long> c_or(n), c_and(n), c_xor(n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                c_or[i | j] += a[i] * b[j], c_and[i & j] += a[i] * b[j], c_xor[i ^ j] += a[i] * b[j];
            }
        }

        assert(fwht::or_convolution(a, b) == c_or);
        assert(fwht::and_convolution(a, b) == c_and);
        assert(fwht::xor_convolution(a, b) == c_xor);
    }

    return 0;
}
