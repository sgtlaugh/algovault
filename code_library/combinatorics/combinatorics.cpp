/***
 *
 * Combinatorics
 * Factorial tables, nCr and nPr modulo a prime
 *
 * Complexity: O(n) to build Combinatorics(n, mod), O(1) per query afterwards
 *
 * Combinatorics comb(n, mod): tables for 0..n, mod must be a prime with n < mod < 2^31
 * comb.inverse(i): inverse of i for 1 <= i <= n, comb.factorial(i), comb.inv_factorial(i)
 * comb.nCr(n, k), comb.nPr(n, k): 0 when k < 0 or k > n, n must not exceed the table size
 *
 * The tables are int to halve their memory, so read them through the calls, which return long long
 * Extended gcd, modular inverse for any modulus and diophantine equations: extended_gcd.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Combinatorics{
    long long mod;
    vector<int> fact, inv_fact, inv;

    Combinatorics(int n, long long mod) : mod(mod), fact(n + 1), inv_fact(n + 1), inv(n + 1){
        assert(n >= 0 && n < mod && mod < (1LL << 31));

        fact[0] = 1;
        for (int i = 1; i <= n; i++) fact[i] = (long long)fact[i - 1] * i % mod;

        if (n >= 1) inv[1] = 1;
        for (int i = 2; i <= n; i++) inv[i] = mod - (mod / i) * inv[mod % i] % mod;  /// needs a prime mod

        inv_fact[0] = 1;
        for (int i = 1; i <= n; i++) inv_fact[i] = (long long)inv_fact[i - 1] * inv[i] % mod;
    }

    long long inverse(int i) const{
        assert(1 <= i && i < (int)inv.size());
        return inv[i];
    }

    long long factorial(int i) const{
        assert(0 <= i && i < (int)fact.size());
        return fact[i];
    }

    long long inv_factorial(int i) const{
        assert(0 <= i && i < (int)fact.size());
        return inv_fact[i];
    }

    long long nCr(int n, int k) const{
        assert(n < (int)fact.size());
        if (k < 0 || k > n) return 0;
        return (long long)fact[n] * inv_fact[k] % mod * inv_fact[n - k] % mod;
    }

    long long nPr(int n, int k) const{
        assert(n < (int)fact.size());
        if (k < 0 || k > n) return 0;
        return (long long)fact[n] * inv_fact[n - k] % mod;
    }
};

int main(){
    Combinatorics comb(10, 1000000007);
    assert(comb.nCr(5, 2) == 10);
    assert(comb.nCr(5, 6) == 0);      /// k > n
    assert(comb.nCr(10, 5) == 252);
    assert(comb.nPr(5, 2) == 20);
    assert(comb.factorial(10) == 3628800);
    assert(comb.inverse(2) == 500000004);
    assert(comb.inverse(10) == 700000005);  /// 10 * 700000005 = 7 * (10^9 + 7) + 1

    Combinatorics small(4, 5);
    assert(small.nCr(4, 2) == 1);         /// 6 mod 5
    assert(small.inv_factorial(4) == 4);  /// 4! = 24 = 4 mod 5, 4 * 4 = 16 = 1 mod 5
    return 0;
}
