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
    if ((g % p) == 0) return (h == 1) ? 0 : -1;

    int i, x, y, z, r, m, counter = 0;
    long long v = 1, d = 1, mul = 1, t = 1 % p;

    for (i = 0; i < 100; i++){
        if (t == h) return i;
        t = (t * g) % p;
    }

    while ((v = __gcd(g, p)) > 1){
        if (h % v) return -1;
        h /= v, p /= v;
        d = (d * (g / v)) % p;
        counter++;
    }

    m = sqrt(p - 0.5) + 1;
    tr1::unordered_map <int, int> mp;

    for (i = 0; i < m; i++){
        if (!mp[mul]) mp[mul] = i + 1;
        mul = (mul * g) % p;
    }

    for (i = 0; i < m; i++){
        z = ext_euclid(d, p, x, y);
        r = ((((long long)x * h) / z) % p + p) % p;
        if (mp.count(r)) return ((i * m) + mp[r] + counter - 1);
        d = (d * mul) % p;
    }

    return -1;
}

int main(){
    printf("%d\n", discrete_log(1, 1, 1));                            /// -1
    printf("%d\n", discrete_log(2, 1, 3));                            /// 0
    printf("%d\n", discrete_log(2, 3, 4));                            /// -1
    printf("%d\n", discrete_log(5, 33, 58));                          /// 9
    printf("%d\n", discrete_log(3589, 58, 97));                       /// -1
    printf("%d\n", discrete_log(3589, 1, 97));                        /// 0
    printf("%d\n", discrete_log(1000000000, 666666667, 1000000007));  /// 942190576

    return 0;
}
