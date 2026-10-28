#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/bell_numbers.cpp"
#undef main

/// Counts set partitions of n items by growing restricted strings: item i joins a used group or opens group k
unsigned long long count_partitions(int i, int n, int k){
    if (i == n) return 1;
    unsigned long long res = 0;
    for (int g = 0; g <= k; g++) res += count_partitions(i + 1, n, max(k, g + 1));
    return res;
}

/// B_0..B_n modulo m as row sums of S(i, k) = k S(i - 1, k) + S(i - 1, k - 1), independent of the Bell triangle
vector<long long> stirling_sums(int n, long long m){
    vector<unsigned __int128> s = {1};
    vector<long long> res = {(long long)(1 % m)};
    for (int i = 1; i <= n; i++){
        s.push_back(0);
        for (int k = i; k > 0; k--) s[k] = ((unsigned __int128)k * s[k] + s[k - 1]) % m;
        s[0] = 0;
        unsigned __int128 sum = 0;
        for (auto x : s) sum = (sum + x) % m;
        res.push_back((long long)sum);
    }
    return res;
}

/// B_0..B_n modulo a prime p by iterating Touchard's recurrence B_(i + p) = B_i + B_(i + 1) from B_0..B_(p - 1)
vector<long long> touchard_iterate(int n, long long p){
    vector<long long> b = stirling_sums(p - 1, p);
    while ((int)b.size() <= n) b.push_back((b[b.size() - p] + b[b.size() - p + 1]) % p);
    b.resize(n + 1);
    return b;
}

int main(){
    /// Exact counts from the brute force, read through moduli larger than every value
    for (int n = 0; n <= 10; n++){
        long long exact = count_partitions(0, n, 0);
        assert(bell_numbers(n, 1LL << 62)[n] == exact);
        assert(bell_numbers_ntt(n)[n] == exact);
        assert(bell_number(n, 1048573) == exact);
    }

    /// Small and extreme moduli, m = 1 and the 2^62 bound included, against the Stirling sums
    vector<long long> moduli = {1, 2, 3, 4, 6, 1000000007, 998244353, (1LL << 62) - 1, 1LL << 62};
    for (long long it = 0; it < stress::scaled(10); it++){
        moduli.push_back(stress::rand_int(1, 1LL << 62));
        moduli.push_back(stress::rand_int(1, 1000));
    }
    for (long long m : moduli){
        int n = stress::rand_int(0, 400);
        assert(bell_numbers(n, m) == stirling_sums(n, m));
    }

    /// The NTT path against the triangle: the one-NTT prime, the CRT path, the largest prime below 2^30, MOD = n + 3
    auto b1 = bell_numbers(2000, 998244353), b2 = bell_numbers(2000, 1000000007), b3 = bell_numbers(2000, 1073741789);
    for (long long it = 0; it < stress::scaled(5); it++){
        int n = it ? stress::rand_int(0, 2000) : 2000;
        assert(bell_numbers_ntt(n) == vector<long long>(b1.begin(), b1.begin() + n + 1));
        assert(bell_numbers_ntt<1000000007>(n) == vector<long long>(b2.begin(), b2.begin() + n + 1));
        assert(bell_numbers_ntt<1073741789>(n) == vector<long long>(b3.begin(), b3.begin() + n + 1));
    }
    assert(bell_numbers_ntt<2003>(2000) == bell_numbers(2000, 2003));
    auto big = bell_numbers_ntt(1 << 14);
    auto head = stirling_sums(1000, 998244353);
    assert(vector<long long>(big.begin(), big.begin() + 1001) == head);

    /// Every n around and past p for each small prime, against the iterated recurrence
    vector<long long> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 97};
    for (long long p : primes){
        auto expected = touchard_iterate(4 * p + 10, p);
        for (int n = 0; n <= 4 * p + 10; n++) assert(bell_number(n, p) == expected[n]);
    }

    /// Random far n against the iterated recurrence, which walks every term up to n
    for (long long it = 0; it < stress::scaled(10); it++){
        long long p = primes[stress::rand_int(0, primes.size() - 1)];
        int n = stress::rand_int(0, 200000);
        assert(bell_number(n, p) == touchard_iterate(n, p)[n]);
    }

    /// n up to 1e18: B_n mod p has period (p^p - 1) / (p - 1), 3, 13 and 781 for p = 2, 3, 5
    vector<pair<long long, long long>> periods = {{2, 3}, {3, 13}, {5, 781}};
    for (auto [p, period] : periods){
        auto expected = touchard_iterate(period, p);
        for (long long it = 0; it < stress::scaled(20); it++){
            long long n = it ? stress::rand_int(0, (long long)1e18) : (long long)1e18;
            assert(bell_number(n, p) == expected[n % period]);
        }
    }

    /// Touchard's congruence itself at huge n for primes too big to iterate up to n
    for (long long p : {101LL, 401LL}){
        long long n = stress::rand_int((long long)1e18 - 1000, (long long)1e18 - p - 1);
        assert(bell_number(n + p, p) == (bell_number(n, p) + bell_number(n + 1, p)) % p);
    }

    return 0;
}
