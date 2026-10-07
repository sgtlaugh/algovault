#include "common.h"

#define main library_main
#include "../code_library/mobius_function.cpp"
#undef main

int main(){
    /// mu from smallest prime factors, independent of both library sieves
    vector<int> spf(MAX, 0);
    vector<signed char> expected(MAX, 0);
    for (int i = 2; i < MAX; i++){
        if (spf[i]) continue;
        for (long long j = i; j < MAX; j += i) if (!spf[j]) spf[j] = i;
    }
    expected[1] = 1;
    for (int i = 2; i < MAX; i++){
        int p = spf[i], rest = i / p;
        expected[i] = rest % p == 0 ? 0 : -expected[rest];
    }

    generate_mobius_fast();
    assert(equal(expected.begin(), expected.end(), mu));
    assert(len == 664579);

    memset(mu, 0, sizeof(mu));
    generate_mobius();
    assert(equal(expected.begin(), expected.end(), mu));
    return 0;
}
