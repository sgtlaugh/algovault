/***
 *
 * Necklaces
 * Counts rings of n beads in k colors up to rotation (necklaces), or up to rotation and turning over (bracelets)
 *
 * Complexity:
 *   - Any number of each color: O(sqrt(n) + d(n) log n) for d(n) divisors of n, at most 6720 for n <= 10^12
 *   - Fixed counts: O(sqrt(g) + d(g) (n + k)) for k colors, g the gcd of the counts
 *
 * necklaces<T>(n, k): number of necklaces, exact in T
 * bracelets<T>(n, k): number of bracelets, exact in T
 * necklaces(n, k, m): number of necklaces modulo m, needs n m < 2^64
 * bracelets(n, k, m): number of bracelets modulo m, needs 2 n m < 2^64
 * necklaces<T>(counts): number of necklaces with exactly counts[j] beads of color j, exact in T
 * bracelets<T>(counts): number of bracelets with exactly counts[j] beads of color j, exact in T
 *
 * Notes:
 *   - n >= 1, k >= 0 and m >= 1, and m may be any modulus, prime or not, even one dividing n
 *   - T is any integer type with +, * and / by a long long: long long, __int128, Bignum from bignum.cpp
 *   - T must hold n times the necklace count and 2 n times the bracelet count
 *   - counts: n beads in all, n >= 1, and T must hold 2 n times the n! / (counts[0]! counts[1]! ...) arrangements
 *   - counts modulo m: reduce the Bignum result, since the multinomials need exact division
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace burnside{
    /// A residue modulo mod = n m, whose division by n is exact on the integers and gives the answer modulo m
    struct Mod{
        unsigned long long v, mod;

        Mod operator + (const Mod& o) const{
            return {v >= mod - o.v ? v - (mod - o.v) : v + o.v, mod};
        }

        Mod operator * (const Mod& o) const{
            return {(unsigned long long)((unsigned __int128)v * o.v % mod), mod};
        }

        Mod operator * (long long x) const{
            return *this * Mod{(unsigned long long)x % mod, mod};
        }

        Mod operator / (long long x) const{
            return {v / x, mod};
        }
    };

    /// Left to right, so no intermediate exceeds b^e and an exact T never overflows before the answer does
    template<typename T>
    T power(const T& b, long long e){
        T res = b;
        for (int bit = __lg(e) - 1; bit >= 0; bit--){
            res = res * res;
            if (e >> bit & 1) res = res * b;
        }

        return res;
    }

    /// Every divisor d of n with phi(d)
    vector<pair<long long, long long>> divisors(long long n){
        vector<pair<long long, long long>> divs = {{1, 1}};
        long long rest = n;
        for (long long p = 2; rest > 1; p++){
            if (p > rest / p) p = rest;
            if (rest % p) continue;

            int len = divs.size();
            for (long long pj = 1; rest % p == 0; pj *= p){
                rest /= p;
                for (int i = 0; i < len; i++){
                    long long d = divs[i].first, phi = divs[i].second;
                    divs.push_back({d * pj * p, phi * pj * (p - 1)});
                }
            }
        }

        return divs;
    }

    /// (sum of c / len)! / prod (c / len)! over the counts c, rounding down, by exact divisions in T
    template<typename T>
    T multinomial(const vector<int>& counts, long long len){
        T res(1);
        long long total = 0;
        for (int c : counts){
            for (long long i = 1; i <= c / len; i++){
                total++;
                res = res * total / i;
            }
        }

        return res;
    }

    /// Sum of phi(d) k^(n / d) over the divisors d of n, n times the necklace count
    template<typename T>
    T rotations(long long n, const T& k){
        vector<pair<long long, long long>> divs = divisors(n);
        T res = power(k, n);
        for (size_t i = 1; i < divs.size(); i++) res = res + power(k, n / divs[i].first) * divs[i].second;
        return res;
    }

    /// Rotations with cycles of length d fix arrangements only when d divides every count, so d runs over g's divisors
    template<typename T>
    T rotations(const vector<int>& counts){
        long long g = 0;
        for (int c : counts){
            assert(c >= 0);
            g = gcd(g, (long long)c);
        }
        assert(g >= 1);

        T res(0);
        for (auto [d, phi] : divisors(g)) res = res + multinomial<T>(counts, d) * phi;
        return res;
    }

    /// Sum of colorings fixed by each reflection, 2 n times the bracelet count minus rotations(n, k)
    template<typename T>
    T reflections(long long n, const T& k){
        if (n % 2) return power(k, (n + 1) / 2) * n;
        T half = power(k, n / 2);
        return (half + half * k) * (n / 2);
    }
}

template<typename T>
T necklaces(long long n, long long k){
    assert(n >= 1 && k >= 0);
    return burnside::rotations(n, T(k)) / n;
}

template<typename T>
T bracelets(long long n, long long k){
    assert(n >= 1 && k >= 0 && n < (1LL << 62));
    T kt(k);
    return (burnside::rotations(n, kt) + burnside::reflections(n, kt)) / (2 * n);
}

long long necklaces(long long n, long long k, long long m){
    assert(n >= 1 && k >= 0 && m >= 1 && (unsigned __int128)n * m <= ULLONG_MAX);
    unsigned long long mod = (unsigned long long)n * m;
    return (burnside::rotations(n, burnside::Mod{k % mod, mod}) / n).v;
}

long long bracelets(long long n, long long k, long long m){
    assert(n >= 1 && k >= 0 && m >= 1 && (unsigned __int128)2 * n * m <= ULLONG_MAX);
    unsigned long long mod = 2ULL * n * m;
    burnside::Mod kt{k % mod, mod};
    return ((burnside::rotations(n, kt) + burnside::reflections(n, kt)) / (2 * n)).v;
}

template<typename T>
T necklaces(const vector<int>& counts){
    T res = burnside::rotations<T>(counts);
    return res / accumulate(counts.begin(), counts.end(), 0LL);
}

/// Every reflection fixes n times multinomial(counts, 2) colorings in all, when at most 2 counts are odd, else none
template<typename T>
T bracelets(const vector<int>& counts){
    T res = burnside::rotations<T>(counts);
    long long n = accumulate(counts.begin(), counts.end(), 0LL), odd = 0;
    for (int c : counts) odd += c % 2;

    if (odd <= 2) res = res + burnside::multinomial<T>(counts, 2) * n;
    return res / (2 * n);
}

int main(){
    /// 4 beads in 2 colors: 0000 0001 0011 0101 0111 1111, and turning over merges no two of them
    assert(necklaces<long long>(4, 2) == 6);
    assert(bracelets<long long>(4, 2) == 6);

    /// 6 beads in 3 colors: 130 necklaces, 38 pairs of them mirror images of each other
    assert(necklaces<long long>(6, 3) == 130);
    assert(bracelets<long long>(6, 3) == 92);

    /// 7 is prime: (2^7 - 2) / 7 + 2 = 20 necklaces, and each reflection fixes 2^4 colorings, (140 + 112) / 14 = 18
    assert(necklaces(7, 2, 1000000007) == 20);
    assert(bracelets(7, 2, 1000000007) == 18);

    /// m = 3 divides n = 6, where dividing by the inverse of n is impossible
    assert(necklaces(6, 3, 3) == 130 % 3);
    assert(bracelets(6, 3, 3) == 92 % 3);

    /// 6 beads, 2 of each of 3 colors: rotating by 3 fixes 3! = 6 of the 90 arrangements, (90 + 6) / 6 = 16
    assert(necklaces<long long>({2, 2, 2}) == 16);
    /// 6 reflections fix 6 each: 3 colors for the 2 beads on the axis times 2, or 3! for 3 pairs, (96 + 36) / 12 = 11
    assert(bracelets<long long>({2, 2, 2}) == 11);

    return 0;
}
