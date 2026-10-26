#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/miller_rabin.cpp"
#undef main

/// Strong probable prime test of odd n > 2 to base a, with plain __int128 arithmetic
bool reference_passes(unsigned long long n, unsigned long long a){
    unsigned long long d = n - 1;
    int s = 0;
    while (!(d & 1)) d >>= 1, s++;
    auto mul = [&](unsigned long long x, unsigned long long y){ return (unsigned long long)((unsigned __int128)x * y % n); };

    unsigned long long x = 1, b = a % n;
    for (unsigned long long e = d; e; e >>= 1, b = mul(b, b)) if (e & 1) x = mul(x, b);
    if (x == 1 || x == n - 1) return true;
    for (int r = 1; r < s; r++) if ((x = mul(x, x)) == n - 1) return true;
    return false;
}

/// Independent reference: Miller-Rabin with the first 12 prime bases is exact below 3.18 * 10^23
bool reference_is_prime(unsigned long long n){
    if (n < 2) return false;
    const int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : bases) if (n % p == 0) return n == (unsigned long long)p;
    for (int a : bases) if (!reference_passes(n, a)) return false;
    return true;
}

/// One witness among liars in a random lane must be rejected, so a lane whose verdict is dropped or misread fails
template<size_t K>
void check_lanes(long long n, const vector<int>& liars, const vector<int>& witnesses){
    array<int, K> bases;
    for (auto& a : bases) a = liars[stress::rand_int(0, liars.size() - 1)];
    assert(prm::miller_rabin(n, bases));

    bases[stress::rand_int(0, K - 1)] = witnesses[stress::rand_int(0, witnesses.size() - 1)];
    assert(!prm::miller_rabin(n, bases));
}

int main(){
    const int N = 200000;
    vector<bool> sieve(N, true);
    sieve[0] = sieve[1] = false;
    for (int i = 2; i * i < N; i++) if (sieve[i]) for (int j = i * i; j < N; j += i) sieve[j] = false;

    for (int n = -5; n < N; n++) assert(prm::is_prime(n) == (n >= 0 && sieve[n]));

    /// Composites passing base 2 with no prime factor below 80, so trial division leaves them to the lanes
    for (int n = 81; n < N; n += 2){
        bool small_factor = false;
        for (int p = 3; p < 80; p += 2) small_factor |= n % p == 0;
        if (sieve[n] || small_factor || !reference_passes(n, 2)) continue;

        /// Bases congruent to 1 or -1 pass every n, and a base above n must be reduced first
        vector<int> liars = {1, n - 1, n + 1, 2 * n - 1}, witnesses;
        for (int a = 2; a < 3000; a++) (reference_passes(n, a) ? liars : witnesses).push_back(a);
        for (int it = 0; it < stress::scaled(300); it++){
            check_lanes<3>(n, liars, witnesses);
            check_lanes<7>(n, liars, witnesses);
        }
    }

    /// Strong pseudoprimes to small base sets, 4759123141 to {2, 7, 61} so is_prime must switch to 7 bases there
    const long long special[] = {2047, 1373653, 25326001, 3215031751LL, 2152302898747LL, 3474749660383LL, 341550071728321LL,
                                 3825123056546413051LL, 561, 1105, 1729, 299210837, 407521, 2147483647, 2147483649LL,
                                 4294967291LL, 4294967295LL, 4294967311LL, 4759123141LL, 9223372036854775783LL, LLONG_MAX};
    for (long long n : special) assert(prm::is_prime(n) == reference_is_prime(n));

    for (long long it = 0; it < stress::scaled(1000000); it++){
        long long n = (long long)(stress::rng()() >> stress::rand_int(1, 62));
        assert(prm::is_prime(n) == reference_is_prime(n));
    }

    /// Semiprimes of two large primes are the hardest composites for Miller-Rabin, trial division cannot reject them
    auto random_prime = [](){
        long long p;
        do p = stress::rand_int(65537, 3000000000LL); while (!reference_is_prime(p));
        return p;
    };

    for (long long it = 0; it < stress::scaled(10000); it++){
        long long p = random_prime(), q = random_prime();
        if ((unsigned __int128)p * q > LLONG_MAX) continue;
        assert(!prm::is_prime(p * q));
    }

    return 0;
}
