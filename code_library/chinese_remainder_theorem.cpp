#include <bits/stdtr1c++.h>

using namespace std;


long long extended_gcd(long long a, long long b, long long& x, long long& y){
    /// Bezout's identity, ax + by = gcd(a,b)

    if (!b){
        y = 0, x = 1;
        return a;
    }

    long long g = extended_gcd(b, a % b, y, x);
    y -= ((a / b) * x);
    return g;
}

long long mod_inverse(long long a, long long m){
    /// inverse exists if and only if a and m are co-prime

    long long x, y, inv;
    extended_gcd(a, m, x, y);
    inv = (x + m) % m;
    return inv;
}

/***
 * finds the unique solution x modulo M (product of mods) for which x % mods[i] = ar[i]
 * mods must be pairwise co-prime
***/

long long CRT(int n, long long* ar, long long* mods){
    long long x, y, res = 0, M = 1;

    for (int i = 0; i < n; i++) M *= mods[i];
    for (int i = 0; i < n; i++){
        x = M / mods[i];
        y = mod_inverse(x, mods[i]);
        res = (res + (((x * ar[i]) % M) * y)) % M;
    }

    return res;
}

int main(){
    long long mods[] = {3, 5, 7};
    long long ar[] = {2, 3, 2};
    long long res = CRT(3, ar, mods);

    printf("%lld\n", res);  /// 23
    return 0;
}
