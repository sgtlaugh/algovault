/***
 *
 * Returns the pisano period of fibonacci sequence modulo n
 * fib(0) = 0, fib(1) = 1, fib(x) = fib(x - 1) + fib(x - 2) for x > 1
 * If each term of the fibonacci sequence is taken modulo m, then the sequence repeats after every k terms
 * The length of this cycle is known as the pisano period
 *
 * For example,
 * First few fibonacci numbers:  0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55
 * First few fibonacci modulo 4: 0, 1, 1, 2, 3, 1, 0, 1,  1,  2,  3
 * the pisano period for m = 4 is therefore 6, as the sequence repeats after every 6 terms
 *
 * Supports n up to 1.5 * 10^18, the period can be as large as 6n
 * Requires __int128 (64-bit GCC or Clang)
 *
 * Don't forget to initialize pollard by by calling rho::init() before using
 *
***/

#include <bits/stdc++.h>

using namespace std;

typedef unsigned long long u64;
typedef unsigned __int128 u128;

/***
 *
 * Montgomery multiplication for a fixed odd modulus n, values are stored as x * 2^64 mod n
 * Needs no division, ~1.6x faster than the long double trick and ~4x faster than __int128 % n
 * The modulus must be odd, even moduli give wrong results silently
 *
***/

struct Montgomery{
    u64 n, inv, r2;

    Montgomery(u64 n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// Newton iteration, doubles the correct bits of n^-1 mod 2^64 each step
        r2 = -n % n;
        r2 = (u128)r2 * r2 % n;
    }

    inline u64 reduce(u128 x) const{
        u64 q = (u64)x * inv, m = ((u128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    inline u64 mul(u64 x, u64 y) const{
        return reduce((u128)x * y);
    }

    inline u64 add(u64 x, u64 y) const{
        return (x += y) >= n ? x - n : x;
    }

    inline u64 to_mont(u64 x) const{
        return mul(x, r2);
    }

    inline u64 from_mont(u64 x) const{
        return reduce(x);
    }

    inline u64 pow(u64 x, u64 e) const{
        u64 res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }
        return res;
    }
};

inline unsigned long long gcd(unsigned long long u, unsigned long long v){
    if (!u || !v) return u | v;
    if (u == 1 || v == 1) return 1;

    int shift = __builtin_ctzll(u | v);
    u >>= __builtin_ctzll(u);
    do{
        v >>= __builtin_ctzll(v);
        if (u > v) swap(u, v);
        v = v - u;
    } while (v);

    return u << shift;
}

inline long long lcm(long long a, long long b){
    return (a / gcd(a, b)) * b;
}

inline void fib(u64& x, u64& y, long long n, const Montgomery& mont){
    if (!n) x = 0, y = mont.to_mont(1);
    else{
        u64 a, b;
        fib(a, b, n >> 1, mont);
        u64 z = mont.add(b, b) + mont.n - a;
        if (z >= mont.n) z -= mont.n;

        x = mont.mul(a, z);
        y = mont.add(mont.mul(a, a), mont.mul(b, b));

        if (n & 1){
            x = mont.add(x, y);
            swap(x, y);
        }
    }
}

/// m must be odd
inline pair<long long, long long> fib(long long n, long long m){
    Montgomery mont(m);
    u64 x, y;
    fib(x, y, n, mont);
    return pair<long long, long long>(mont.from_mont(x), mont.from_mont(y));
}

namespace rho{
    const int MAXP = 1000010;
    const int BASE[] = {2, 450775, 1795265022, 9780504, 28178, 9375, 325};

    int primes[MAXP], spf[MAXP];
    long long divisors[130172];

    inline bool miller_rabin(long long n){
        if (n <= 2 || !(n & 1)) return n == 2;
        if (n < MAXP) return spf[n] == n;

        long long s = 0, r = n - 1;
        for (; !(r & 1); r >>= 1, s++) {}

        Montgomery mont(n);
        u64 c, d, one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        for (int i = 0; i < 7; i++){
            long long a = BASE[i] % n;
            if (!a) continue;  /// n divides the base (e.g. 299210837), this base says nothing about n
            c = mont.pow(mont.to_mont(a), r);
            for (int j = 0; j < s; j++){
                d = mont.mul(c, c);
                if (d == one && c != one && c != minus_one) return false;
                c = d;
            }

            if (c != one) return false;
        }
        return true;
    }

    inline void init(){
        int i, j, k, cnt = 0;

        for (i = 2; i < MAXP; i++){
            if (!spf[i]) primes[cnt++] = spf[i] = i;
            for (j = 0; (k = i * primes[j]) < MAXP; j++){
                spf[k] = primes[j];
                if(spf[i] == spf[k]) break;
            }
        }
    }

    /// Brent's cycle detection, one gcd per BATCH steps instead of one per step
    long long pollard_rho(long long n){
        if (!(n & 1)) return 2;

        const int BATCH = 128;
        Montgomery mont(n);
        auto next = [&](u64 x, u64 c){ return mont.add(mont.mul(x, x), c); };
        auto dist = [](u64 x, u64 y){ return x > y ? x - y : y - x; };

        while (1){
            u64 c = rand() % (n - 1) + 1, y = rand() % n, x = y, ys = y, q = mont.to_mont(1), g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y, c);

                for (long long k = 0; k < r && g == 1; k += BATCH){
                    ys = y;
                    for (long long i = 0; i < BATCH && i < r - k; i++){
                        y = next(y, c);
                        q = mont.mul(q, dist(x, y));
                    }
                    g = gcd(q, n);
                }
            }

            /// The batch overshot to a multiple of n, replay it one step at a time
            if (g == (u64)n){
                do{
                    ys = next(ys, c);
                    g = gcd(dist(x, ys), n);
                } while (g == 1);
            }

            if (g != (u64)n) return g;
        }
    }

    vector <long long> factorize(long long n){
        if (n == 1) return vector <long long>();
        if (miller_rabin(n)) return vector<long long> {n};

        vector <long long> v, w;
        while (n > 1 && n < MAXP){
            v.push_back(spf[n]);
            n /= spf[n];
        }

        if (n >= MAXP) {
            long long x = pollard_rho(n);
            v = factorize(x);
            w = factorize(n / x);
            v.insert(v.end(), w.begin(), w.end());
        }

        sort(v.begin(), v.end());
        return v;
    }

    vector <long long> get_divisors(long long n){
        int j, k, l, len = 0;

        auto factors = factorize(n);
        map <long long, int> prime_count;
        for (auto x: factors) prime_count[x]++;

        divisors[len++] = 1;
        for (auto it: prime_count){
            long long x = it.first;
            for (k = len, j = 0; j < it.second; j++, x *= it.first){
                for (l = 0; l < k; l++){
                    divisors[len++] = x * divisors[l];
                }
            }
        }

        sort(divisors, divisors + len);
        return vector <long long>(divisors, divisors + len);
    }
}

long long pi(long long p){
    if (p == 2) return 3;
    if (p == 5) return 20;

    int flag = p % 10;
    long long n = (flag == 1 || flag == 9) ? p - 1 : 2 * (p + 1);

    auto divisors = rho::get_divisors(n);
    for (auto d: divisors){
        auto seq = fib(d, p);
        if (seq.first == 0 && seq.second == 1) return d;
    }

    return -1;
}

long long pisano_period(long long n){
    auto pfs = rho::factorize(n);
    unordered_map <long long, int> pf_cnts;
    for (auto pf: pfs) pf_cnts[pf]++;

    long long res = 1;
    for (auto it: pf_cnts){
        long long x = 1;
        for (int i = 1; i < it.second; i++) x *= it.first;
        res = lcm(res, pi(it.first) * x);
    }

    return res;
}

int main(){
    rho::init();

    assert(pisano_period(1) == 1);
    assert(pisano_period(4) == 6);
    assert(pisano_period(47) == 32);
    assert(pisano_period(1260) == 240);
    assert(pisano_period(510510) == 5040);
    assert(pisano_period(2147483647LL) == 4294967296LL);
    assert(pisano_period(2147483648LL) == 3221225472LL);
    assert(pisano_period(133049351085651000LL) == 5212409202000LL);
    assert(pisano_period(1000000000000000000LL) == 1500000000000000000LL);

    return 0;
}
