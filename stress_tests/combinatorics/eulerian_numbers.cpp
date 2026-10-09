#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/eulerian_numbers.cpp"
#undef main

long long power(long long x, long long n, long long p){
    long long res = 1;
    for (x %= p; n; n >>= 1, x = x * x % p) if (n & 1) res = res * x % p;
    return res;
}

/// A(n, k) = sum (-1)^j C(n + 1, j) (k + 1 - j)^n over j <= k, checked for the whole row n modulo a prime p
void check_row(const EulerianNumbers& eul, int n, long long p){
    vector<long long> binom(n + 2), pw(n + 2);
    binom[0] = 1;
    for (int j = 1; j <= n + 1; j++) binom[j] = binom[j - 1] * (n + 2 - j) % p * power(j, p - 2, p) % p;
    for (int x = 0; x <= n + 1; x++) pw[x] = power(x, n, p);

    for (int k = 0; k < n; k++){
        long long a = 0;
        for (int j = 0; j <= k; j++) a = (a + (j & 1 ? p - 1 : 1) * binom[j] % p * pw[k + 1 - j]) % p;
        assert(eul.get(n, k) == a);
    }
    assert(eul.get(n, n) == (n == 0));
}

int main(){
    const int MAXN = 5009;
    vector<int> moduli = {1, 2, 7, 1000000007, 2147483647};
    for (long long it = 0; it < stress::scaled(10); it++) moduli.push_back(stress::rand_int(1, 2147483647));

    /// Ascents counted over every permutation
    for (int m : moduli){
        EulerianNumbers eul(8, m);
        for (int n = 1; n <= 8; n++){
            vector<int> p(n), count(n + 1, 0);
            iota(p.begin(), p.end(), 1);
            do {
                int ascents = 0;
                for (int i = 1; i < n; i++) ascents += p[i] > p[i - 1];
                count[ascents]++;
            } while (next_permutation(p.begin(), p.end()));
            for (int k = 0; k <= n; k++) assert(eul.get(n, k) == count[k] % m);
        }
        assert(eul.get(0, 0) == 1 % m);
    }

    for (long long p : {1000000007LL, 2147483647LL}){
        EulerianNumbers eul(MAXN, p);
        for (int n = 0; n <= 60; n++) check_row(eul, n, p);
        check_row(eul, MAXN, p);
        for (long long it = 0; it < stress::scaled(3); it++) check_row(eul, stress::rand_int(it % 2 ? 61 : MAXN - 100, MAXN), p);
    }

    return 0;
}
