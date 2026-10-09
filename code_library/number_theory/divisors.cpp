/***
 *
 * Divisors of a single number
 *
 * Complexity: O(sqrt(n)) to factorize by trial division, O(d(n) log d(n)) to list and sort the divisors
 *
 * divisors(n): every divisor of 1 <= n <= 1e14 in ascending order
 * divisors_from_factors(f): the same from prime factors with multiplicity, e.g. pollard_rho.cpp's factorize(n)
 *     use this for n up to 2^63, products stay below n so nothing overflows
 *
 * Tables of the divisor count or divisor sum for every n up to a limit: linear_sieve.cpp's multiplicative()
 * Every divisor of every n up to a limit: all_divisors.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> divisors_from_factors(vector<long long> factors){
    sort(factors.begin(), factors.end());
    vector<long long> res = {1};

    for (size_t i = 0; i < factors.size(); ){
        size_t j = i;
        while (j < factors.size() && factors[j] == factors[i]) j++;
        size_t old = res.size();
        long long power = 1;
        for (size_t k = i; k < j; k++){
            power *= factors[i];
            for (size_t t = 0; t < old; t++) res.push_back(res[t] * power);
        }
        i = j;
    }

    sort(res.begin(), res.end());
    return res;
}

vector<long long> divisors(long long n){
    assert(1 <= n && n <= 100000000000000LL);
    vector<long long> factors;

    for (long long p = 2; p * p <= n; p++){
        for (; n % p == 0; n /= p) factors.push_back(p);
    }
    if (n > 1) factors.push_back(n);

    return divisors_from_factors(factors);
}

int main(){
    assert((divisors(1) == vector<long long>{1}));
    assert((divisors(12) == vector<long long>{1, 2, 3, 4, 6, 12}));
    assert((divisors(97) == vector<long long>{1, 97}));
    assert((divisors(36) == vector<long long>{1, 2, 3, 4, 6, 9, 12, 18, 36}));
    assert(divisors(720720).size() == 240);

    assert((divisors_from_factors({}) == vector<long long>{1}));
    assert((divisors_from_factors({1000000007, 998244353}) == vector<long long>{1, 998244353, 1000000007, 998244353LL * 1000000007}));
    assert(divisors_from_factors(vector<long long>(62, 2)).back() == (1LL << 62));

    return 0;
}
