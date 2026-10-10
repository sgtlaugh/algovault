// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/enumerate_primes
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/fast_sieve.cpp"
#undef main

void write_int(uint32_t x, char end){
    char buf[12];
    int len = 0;
    do buf[len++] = '0' + x % 10; while (x /= 10);
    while (len) putchar(buf[--len]);
    putchar(end);
}

int main(){
    uint32_t n, a, b;
    if (scanf("%u %u %u", &n, &a, &b) != 3) return 0;

    fast_sieve();
    uint32_t pi = upper_bound(primes, primes + prime_cnt, n) - primes;
    uint32_t x = b < pi ? (pi - 1 - b) / a + 1 : 0;

    printf("%u %u\n", pi, x);
    for (uint32_t i = 0; i < x; i++) write_int(primes[a * i + b], i + 1 == x ? '\n' : ' ');
    if (!x) putchar('\n');
    return 0;
}
