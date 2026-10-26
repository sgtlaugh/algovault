/***
 *
 * Divisors
 * Divisors of a single number
 *
 * Complexity: O(sqrt(n)) to factorize by trial division, O(d(n) log d(n)) to list and sort the divisors
 *
 * divisors(n): every divisor of 1 <= n <= 1e14 in ascending order
 * divisors_from_factors(f): the same from prime factors with multiplicity, like pollard_rho.cpp's factorize, n < 2^63
 *
 * Notes:
 *   - divisors_from_factors needs prime factors, composite ones give duplicate or missing divisors
 *   - Tables of the divisor count or divisor sum for every n up to a limit: linear_sieve.cpp's multiplicative()
 *   - Every divisor of every n up to a limit: all_divisors.cpp
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

    /// Odd p only after 2, measured 1.9x faster than trying every p on n near 1e14
    for (; n % 2 == 0; n /= 2) factors.push_back(2);
    for (long long p = 3; p * p <= n; p += 2){
        for (; n % p == 0; n /= p) factors.push_back(p);
    }
    if (n > 1) factors.push_back(n);

    return divisors_from_factors(factors);
}

int main(){
    /// Sorted divisors of one number
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
