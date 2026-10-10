// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/counting_squarefrees
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/du_sieve.cpp"
#undef main
#define main integer_root_main
#include "../../code_library/number_theory/integer_root.cpp"
#undef main

int main(){
    long long n;
    if (scanf("%lld", &n) != 1) return 0;

    /// sum of mu(d) * floor(n / d^2) over d <= sqrt(n), d <= n^(2/5) one by one, the rest in blocks of equal floor(n / d^2)
    long long s = isqrt(n), direct = min(s, max(1LL, (long long)pow((double)n, 0.4)));
    DuSieve du(s, direct);

    long long res = 0;
    for (long long d = 1; d <= direct; d++) res += (du.mertens(d) - du.mertens(d - 1)) * (n / (d * d));

    long long prev = du.mertens(direct);
    for (long long d = direct + 1; d <= s;){
        long long k = n / (d * d), hi = min(s, (long long)isqrt(n / k));
        long long cur = du.mertens(hi);
        res += k * (cur - prev);
        prev = cur, d = hi + 1;
    }

    printf("%lld\n", res);
    return 0;
}
