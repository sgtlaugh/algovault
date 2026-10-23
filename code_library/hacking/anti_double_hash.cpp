/***
 *
 * Anti Double Hash
 * Builds two different strings with the same double polynomial hash, to hack solutions or break weak tests
 *
 * Complexity:
 *   - AntiHash: 2^23 modular multiplies to fill a 4 MB table, independent of the input
 *   - solve: expected O(sqrt(mod1 mod2)) calls to an O(1) hash, the tail plus cycle length of its rho
 *
 * AntiHash hack(base1, mod1, base2, mod2): the attack on h = (h * base + s[i]) % mod under both (base, mod) pairs
 * hack.solve(): {s, t}, two different 64-character strings of '0' and '1' whose hashes match under both pairs
 *
 * Notes:
 *   - solve never terminates if its fixed start value lies on the cycle of the rho
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct AntiHash{
    unsigned long long b[2], m[2], expo[2][64];
    vector<array<array<unsigned long long, 2>, 4>> dp;  /// 4 MB, too big for the stack as a member array

    AntiHash(long long base1, long long mod1, long long base2, long long mod2){
        int l, i, j, k;
        b[0] = base1, m[0] = mod1, b[1] = base2, m[1] = mod2;

        for (l = 0; l < 2; l++){
            for (expo[l][0] = 1, i = 1; i < 64; i++) expo[l][i] = multiply(expo[l][i - 1], b[l], m[l]);
        }

        dp.assign(1 << 16, {});
        for (j = 0; j < (1 << 16); j++){
            for (i = 0; i < 4; i++){
                for (l = 0; l < 2; l++){
                    for (k = 0; k < 16; k++){
                        if (j & (1 << k)) dp[j][i][l] = (dp[j][i][l] + multiply(expo[l][16 * i + k], 1056, m[l])) % m[l];
                        else dp[j][i][l] = (dp[j][i][l] + multiply(expo[l][16 * i + k], 1055, m[l])) % m[l];
                    }
                }
            }
        }
    }

    inline long long multiply(long long a, long long b, long long m){
        a %= m, b %= m;
        if (m < 3037000500LL) return a * b % m;

        #ifdef __SIZEOF_INT128__
            return __int128(a) * b % m;
        #endif

        long double x = a;
        x *= b;
        long long c = (long long)(x / m);

        a *= b;
        a -= c * m;
        if (a >= m) a -= m;
        if (a < 0) a += m;
        return a;
    }

    inline unsigned long long get_hash(unsigned long long x){
        unsigned long long i = 0, h1 = 0, h2 = 0;
        while (x){
            h1 += dp[x & 65535][i][0];
            h2 += dp[x & 65535][i][1];
            i++, x >>= 16;
            if (h1 >= m[0]) h1 -= m[0];
            if (h2 >= m[1]) h2 -= m[1];
        }

        return (h1 << 32) ^ h2 ^ 0x5555555555555555LL;
    }

    string to_string(unsigned long long x){
        char str[65];
        for (int i = 63; i >= 0; i--){
            str[63 - i] = (x & (1LL << i)) ? '1' : '0';
        }
        str[64] = 0;

        return string(str);
    }

    inline pair<string, string> solve(){
        unsigned long long u, v, x, y, s;
        x = y = s = 6666666666666666677LL;

        /// Floyd's cycle finding
        while (1){
            x = get_hash(x);
            y = get_hash(get_hash(y));
            if (x == y) break;
        }

        u = get_hash(x), v = get_hash(s);
        while (!(x != s && u == v)){
            x = u, s = v;
            u = get_hash(x), v = get_hash(s);
        }

        assert(x != s && get_hash(x) == get_hash(s));
        return {to_string(x), to_string(s)};
    }
};

int main(){
    clock_t start = clock();

    long long b1 = 1000003;
    long long m1 = 1000000007;
    long long b2 = 1000007;
    long long m2 = 1000000009;

    string s1, s2;
    tie(s1, s2) = AntiHash(b1, m1, b2, m2).solve();

    assert(s1 == "0100000100000100101100000011010101010010010111111110100001010110");
    assert(s2 == "0110010001010011011100101001101001110001010000001001101111101110");

    long long s1h1 = 0, s1h2 = 0, s2h1 = 0, s2h2 = 0;
    for (int i = 0; i < (int)s1.length(); i++){
        s1h1 = (s1h1 * b1 + s1[i]) % m1;
        s1h2 = (s1h2 * b2 + s1[i]) % m2;
    }
    for (int i = 0; i < (int)s2.length(); i++){
        s2h1 = (s2h1 * b1 + s2[i]) % m1;
        s2h2 = (s2h2 * b2 + s2[i]) % m2;
    }

    assert(s1h1 == s2h1 && s1h2 == s2h2);

    puts(s1.c_str());
    puts(s2.c_str());

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));

    return 0;
}
