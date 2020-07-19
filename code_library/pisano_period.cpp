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
 * Don't forget to initialize pollard by by calling rho::init() before using
 *
***/

#include <bits/stdc++.h>

using namespace std;

inline long long mod_add(long long x, long long y, long long m){
    return (x += y) < m ? x : x - m;
}

inline long long mod_mul(long long x, long long y, long long m){
    #if !defined(STRICT_ANSI) && defined(_GLIBCXX_USE_INT128)
        return (__int128)x * y % m;
    #endif

    long long res = x * y - (long long)((long double)x * y / m + 0.5) * m;
    return res < 0 ? res + m : res;
}

inline long long mod_pow(long long x, long long n, long long m){
    long long res = 1 % m;
    for (; n; n >>= 1){
        if (n & 1) res = mod_mul(res, x, m);
        x = mod_mul(x, x, m);
    }

    return res;
}

inline unsigned long long gcd(unsigned long long u, unsigned long long v){
    if (!u || !v) return u | v;
    if (u == 1 || v == 1) return 1;

    int shift = __builtin_ctzll(u | v);
    u >>= __builtin_ctzll(u);
    do{
        v >>= __builtin_ctzll(v);
        if (u > v) v ^= u ^= v ^= u;
        v = v - u;
    } while (v);

    return u << shift;
}

inline long long lcm(long long a, long long b){
    return (a / gcd(a, b)) * b;
}

inline void fib(long long& x, long long& y, long long n, long long m){
    if (!n) x = 0, y = 1;
    else{
        long long a, b;
        fib(a, b, n >> 1, m);
        long long z = (b << 1) - a;
        if (z < 0) z += m;

        x = mod_mul(a, z, m);
        y = mod_add(mod_mul(a, a, m), mod_mul(b, b, m), m);

        if (n & 1){
            x = mod_add(x, y, m);
            x ^= y ^= x ^= y;
        }
    }
}

inline pair<long long, long long> fib(long long n, long long m){
    long long x = 0, y = 1;
    fib(x, y, n, m);
    return pair<long long, long long>(x, y);
}

namespace rho{
    const int MAXP = 1000010;
    const int BASE[] = {2, 450775, 1795265022, 9780504, 28178, 9375, 325};

    int primes[MAXP], spf[MAXP];
    long long seq[MAXP], divisors[130172];

    inline bool miller_rabin(long long n){
        if (n <= 2 || !(n & 1)) return n == 2;
        if (n < MAXP) return spf[n] == n;

        long long c, d, s = 0, r = n - 1;
        for (; !(r & 1); r >>= 1, s++) {}

        for (int i = 0; i < 7; i++){
            c = mod_pow(BASE[i], r, n);
            for (int j = 0; j < s; j++){
                d = mod_mul(c, c, n);
                if (d == 1 && c != 1 && c != (n - 1)) return false;
                c = d;
            }

            if (c != 1) return false;
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

    long long pollard_rho(long long n){
        while (1){
            long long x = rand() % n, y = x, c = rand() % n, u = 1, v, t = 0;
            long long *px = seq, *py = seq;

            while (1){
                *py++ = y = mod_add(mod_mul(y, y, n), c, n);
                *py++ = y = mod_add(mod_mul(y, y, n), c, n);
                if((x = *px++) == y) break;

                v = u;
                u = mod_mul(u, abs(y - x), n);

                if (!u) return gcd(v, n);
                if (++t == 32){
                    t = 0;
                    if ((u = gcd(u, n)) > 1 && u < n) return u;
                }
            }

            if (t && (u = gcd(u, n)) > 1 && u < n) return u;
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
