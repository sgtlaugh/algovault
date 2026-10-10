/***
 *
 * Pisano Period
 * Period of the Fibonacci sequence modulo n, for 1 <= n <= 1.5 * 10^18
 *
 * Complexity: expected O(n^(1/4)) per Pollard rho split, O(log^2 n) Fibonacci doubling steps per prime factor of n,
 *             O(1e6) once to build the sieve
 *
 * fib(0) = 0, fib(1) = 1, fib(x) = fib(x - 1) + fib(x - 2) for x > 1
 * If each term of the fibonacci sequence is taken modulo m, then the sequence repeats after every k terms
 * The length of this cycle is known as the pisano period
 *
 * For example,
 * First few fibonacci numbers:  0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55
 * First few fibonacci modulo 4: 0, 1, 1, 2, 3, 1, 0, 1,  1,  2,  3
 * the pisano period for m = 4 is therefore 6, as the sequence repeats after every 6 terms
 *
 * PisanoPeriod pisano;    sieves smallest prime factors up to 1e6, keep a single instance around
 * pisano.period(n):       the pisano period of n, can be as large as 6n
 * fib(n, m):              (fib(n) mod m, fib(n + 1) mod m) for an odd modulus m
 *
 * Requires __int128 (64-bit GCC or Clang)
 * Embeds checked copies of pollard_rho.cpp's Montgomery and PollardRho to stay standalone
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n
/// BEGIN COPY montgomery from code_library/number_theory/pollard_rho.cpp
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// Newton iteration, doubles the correct bits of n^-1 mod 2^64 each step
        r2 = -n % n;
        r2 = (unsigned __int128)r2 * r2 % n;
    }

    unsigned long long reduce(unsigned __int128 x) const{
        unsigned long long q = (unsigned long long)x * inv, m = ((unsigned __int128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    unsigned long long mul(unsigned long long x, unsigned long long y) const{
        return reduce((unsigned __int128)x * y);
    }

    unsigned long long add(unsigned long long x, unsigned long long y) const{
        return (x += y) >= n ? x - n : x;
    }

    unsigned long long to_mont(unsigned long long x) const{
        return mul(x, r2);
    }

    unsigned long long pow(unsigned long long x, unsigned long long e) const{
        unsigned long long res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }
        return res;
    }
};
/// END COPY montgomery

/// BEGIN COPY pollard_rho from code_library/number_theory/pollard_rho.cpp
struct PollardRho{
    static constexpr int LIMIT = 1000001;
    vector<int> spf, primes;
    mt19937_64 rng;

    PollardRho() : spf(LIMIT, 0), rng(20260105){
        for (int i = 2; i < LIMIT; i++){
            if (!spf[i]) spf[i] = i, primes.push_back(i);
            for (int p : primes){
                if (p > spf[i] || (long long)i * p >= LIMIT) break;
                spf[i * p] = p;
            }
        }
    }

    bool is_prime(long long n) const{
        if (n < LIMIT) return n >= 2 && spf[n] == n;
        if (!(n & 1)) return false;

        static const long long BASES[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
        long long s = 0, r = n - 1;
        for (; !(r & 1); r >>= 1, s++) {}

        Montgomery mont(n);
        unsigned long long one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        for (long long base : BASES){
            long long a = base % n;
            if (!a) continue;  /// n divides the base, this base says nothing about n
            unsigned long long c = mont.pow(mont.to_mont(a), r);
            for (int j = 0; j < s; j++){
                unsigned long long d = mont.mul(c, c);
                if (d == one && c != one && c != minus_one) return false;
                c = d;
            }

            if (c != one) return false;
        }
        return true;
    }

    vector<long long> factorize(long long n){
        assert(n >= 1);
        vector<long long> res;
        collect(n, res);
        sort(res.begin(), res.end());
        return res;
    }

    void collect(long long n, vector<long long>& out){
        if (n < LIMIT){
            for (; n > 1; n /= spf[n]) out.push_back(spf[n]);
            return;
        }

        if (!(n & 1)){
            out.push_back(2);
            return collect(n >> 1, out);
        }

        if (is_prime(n)){
            out.push_back(n);
            return;
        }

        long long d = split(n);
        collect(d, out);
        collect(n / d, out);
    }

    /// A nontrivial factor of an odd composite n
    long long split(long long n){
        const int BATCH = 128;
        Montgomery mont(n);
        auto next = [&](unsigned long long x, unsigned long long c){ return mont.add(mont.mul(x, x), c); };
        auto dist = [](unsigned long long x, unsigned long long y){ return x > y ? x - y : y - x; };

        while (true){
            unsigned long long c = rng() % (n - 1) + 1, y = rng() % n, x = y, ys = y, q = mont.to_mont(1), g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y, c);

                for (long long k = 0; k < r && g == 1; k += BATCH){
                    ys = y;
                    for (long long i = 0; i < BATCH && i < r - k; i++){
                        y = next(y, c);
                        q = mont.mul(q, dist(x, y));
                    }
                    g = __gcd(q, (unsigned long long)n);
                }
            }

            if (g == (unsigned long long)n){  /// the batch overshot to a multiple of n, replay it one step at a time
                do{
                    ys = next(ys, c);
                    g = __gcd(dist(x, ys), (unsigned long long)n);
                } while (g == 1);
            }

            if (g != (unsigned long long)n) return g;
        }
    }
};
/// END COPY pollard_rho

inline void fib(unsigned long long& x, unsigned long long& y, long long n, const Montgomery& mont){
    if (!n) x = 0, y = mont.to_mont(1);
    else{
        unsigned long long a, b;
        fib(a, b, n >> 1, mont);
        unsigned long long z = mont.add(b, b) + mont.n - a;
        if (z >= mont.n) z -= mont.n;

        x = mont.mul(a, z);
        y = mont.add(mont.mul(a, a), mont.mul(b, b));

        if (n & 1){
            x = mont.add(x, y);
            swap(x, y);
        }
    }
}

inline pair<long long, long long> fib(long long n, long long m){
    assert(m & 1);  /// Montgomery form needs an odd modulus
    Montgomery mont(m);
    unsigned long long x, y;
    fib(x, y, n, mont);
    return pair<long long, long long>(mont.reduce(x), mont.reduce(y));
}

struct PisanoPeriod{
    PollardRho rho;

    long long period(long long n){
        auto pfs = rho.factorize(n);
        unordered_map<long long, int> pf_cnts;
        for (auto pf: pfs) pf_cnts[pf]++;

        long long res = 1;
        for (auto it: pf_cnts){
            long long x = 1;
            for (int i = 1; i < it.second; i++) x *= it.first;
            res = lcm(res, prime_period(it.first) * x);
        }

        return res;
    }

    long long prime_period(long long p){
        if (p == 2) return 3;
        if (p == 5) return 20;

        int flag = p % 10;
        long long n = (flag == 1 || flag == 9) ? p - 1 : 2 * (p + 1);

        /// pi(p) divides n, strip one prime factor at a time (with multiplicity) while the rest is still a period
        long long d = n;
        for (auto q: rho.factorize(n)){
            if (fib(d / q, p) == make_pair(0LL, 1LL)) d /= q;
        }
        return d;
    }
};

int main(){
    PisanoPeriod pisano;

    assert(pisano.period(1) == 1);
    assert(pisano.period(2) == 3);
    assert(pisano.period(3) == 8);
    assert(pisano.period(4) == 6);
    assert(pisano.period(5) == 20);
    assert(pisano.period(10) == 60);
    assert(pisano.period(11) == 10);
    assert(pisano.period(47) == 32);
    assert(pisano.period(1260) == 240);
    assert(pisano.period(510510) == 5040);
    assert(pisano.period(2147483647LL) == 4294967296LL);
    assert(pisano.period(2147483648LL) == 3221225472LL);
    assert(pisano.period(133049351085651000LL) == 5212409202000LL);
    assert(pisano.period(1000000000000000000LL) == 1500000000000000000LL);

    assert(fib(10, 1000000007) == make_pair(55LL, 89LL));
    assert(fib(0, 1) == make_pair(0LL, 0LL));

    PisanoPeriod second;
    assert(second.period(1000000000000000000LL) == 1500000000000000000LL && second.period(47) == 32);

    return 0;
}
