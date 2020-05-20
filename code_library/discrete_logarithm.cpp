#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

int ext_euclid(int a, int b, int& x, int& y){
    if (!b){
        y = 0, x = 1;
        return a;
    }

    int g = ext_euclid(b, a % b, y, x);
    y -= ((a / b) * x);
    return g;
}

/***
 *
 * returns smallest x such that (g^x) % p = h
 * function returns -1 if no solution exists
 * otherwise returns x, the discrete log of h with respect to g modulo p
 *
***/

int discrete_log(int g, int h, int p){
    if (h >= p) return -1;
    if (h == 1 || p == 1) return 0;

    int i, x, y, z, m, c;
    long long v, d, t, mul;

    for (i = 1, t = g; (1LL << i) <= p; i++, t = t * g % p){
        if ((t % p) == h) return i;
    }

    for (c = 0, d = 1, v = 1; (v = __gcd(g, p)) > 1; c++){
        h /= v, p /= v;
        d = d * (g / v) % p;
    }

    m = sqrt(p - 0.666) + 1;
    tr1::unordered_map <int, int> mp;
    for (i = 0, mul = 1; i < m; i++, mul = mul * g % p){
        if (!mp.count(mul)) mp[mul] = i + 1;
    }

    for (i = 0; i < m; i++, d = d * mul % p){
        z = ext_euclid(d, p, x, y);
        z = ((((long long)x * h) / z) % p + p) % p;
        if (mp.count(z)) return i * m + mp[z] + c - 1;
    }

    return -1;
}

int main(){
    printf("%d\n", discrete_log(1, 1, 1));                            /// -1
    printf("%d\n", discrete_log(2, 1, 3));                            ///  0
    printf("%d\n", discrete_log(2, 3, 4));                            /// -1
    printf("%d\n", discrete_log(2, 0, 4));                            ///  2
    printf("%d\n", discrete_log(6, 0, 8));                            ///  3
    printf("%d\n", discrete_log(5, 33, 58));                          ///  9
    printf("%d\n", discrete_log(3589, 58, 97));                       /// -1
    printf("%d\n", discrete_log(3589, 1, 97));                        ///  0
    printf("%d\n", discrete_log(1000000000, 666666667, 1000000007));  ///  942190576

    return 0;
}
