// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/enumerate_primes
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/segmented_sieve.cpp"
#undef main

void write_int(long long x, char end){
    char buf[24];
    int len = 0;
    do buf[len++] = '0' + x % 10; while (x /= 10);
    while (len) putchar(buf[--len]);
    putchar(end);
}

int main(){
    long long n, a, b;
    if (scanf("%lld %lld %lld", &n, &a, &b) != 3) return 0;

    const long long WIDTH = 10000000;
    SegmentedSieve sieve(n);
    long long pi = 0;
    vector<long long> picked;

    for (long long lo = 0; lo <= n; lo += WIDTH){
        long long hi = min(n, lo + WIDTH - 1);
        auto window = sieve.window(lo, hi);
        for (long long i = 0; i <= hi - lo; i++){
            if (!window[i]) continue;
            if (pi >= b && (pi - b) % a == 0) picked.push_back(lo + i);
            pi++;
        }
    }

    printf("%lld %d\n", pi, (int)picked.size());
    for (size_t i = 0; i < picked.size(); i++) write_int(picked[i], i + 1 == picked.size() ? '\n' : ' ');
    if (picked.empty()) putchar('\n');
    return 0;
}
