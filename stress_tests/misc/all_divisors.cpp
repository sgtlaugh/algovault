#include "../common.h"

#define main library_main
#include "../../code_library/misc/all_divisors.cpp"
#undef main

vector<int> trial_division(int x){
    vector<int> res;
    for (int d = 1; (long long)d * d <= x; d++){
        if (x % d) continue;
        res.push_back(d);
        if (d != x / d) res.push_back(x / d);
    }
    sort(res.begin(), res.end());
    return res;
}

void check(int x){
    vector<int> got = divisors[x];
    sort(got.begin(), got.end());
    assert(got == trial_division(x));
}

int main(){
    gen();

    for (int x = 1; x <= 100000; x++) check(x);
    for (int x = MAX - 1000; x < MAX; x++) check(x);
    for (int x : {8648640, 9729720, 1 << 23, 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3, 9999991}) check(x);  /// divisor-rich, prime powers, a prime
    for (long long it = 0; it < stress::scaled(20000); it++) check(stress::rand_int(1, MAX - 1));
    return 0;
}
