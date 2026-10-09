/***
 *
 * Stirling Numbers
 * A whole row of Stirling numbers of the first or second kind modulo m, using FFT multiplication
 *
 * Complexity: O(n log^2 n) for stirling_first, O(n log n) for stirling_second
 *
 * stirling_first(n, m): c(n, k) for k = 0..n, the unsigned numbers counting permutations of n with k cycles
 *     coefficients of x (x + 1) ... (x + n - 1), any modulus 1 <= m <= 2^30
 *     the signed numbers are s(n, k) = (-1)^(n - k) c(n, k)
 * stirling_second(n, p): S(n, k) for k = 0..n, the ways to split n labelled items into k non-empty groups
 *     S(n, k) = sum over i of (-1)^i / i! * (k - i)^n / (k - i)!, needs a prime p with n < p <= 2^30
 *
 * The FFT splits values into 15-bit halves, so products stay exact for lengths up to about 2^20
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace stirling_fft{
    typedef complex<long double> cd;

    /// rt[len / 2 + k] is the k-th root for a level of length len, each computed directly so there is no drift
    vector<cd> rt(2, 1);

    void transform(vector<cd>& a, bool invert){
        int n = a.size();
        for (int len = rt.size(); len < n; len <<= 1){
            rt.resize(2 * len);
            for (int k = 0; k < len; k++) rt[len + k] = polar(1.0L, acosl(-1.0L) * k / len);
        }
        for (int i = 1, j = 0; i < n; i++){
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1){
            for (int i = 0; i < n; i += len){
                for (int k = 0; k < len / 2; k++){
                    cd w = invert ? conj(rt[len / 2 + k]) : rt[len / 2 + k];
                    cd u = a[i + k], v = a[i + k + len / 2] * w;
                    a[i + k] = u + v, a[i + k + len / 2] = u - v;
                }
            }
        }
        if (invert){
            for (auto& x : a) x /= n;
        }
    }

    vector<long long> mod_multiply(const vector<long long>& a, const vector<long long>& b, long long m){
        if (a.empty() || b.empty()) return {};
        int need = a.size() + b.size() - 1, n = 1;
        while (n < need) n <<= 1;

        vector<cd> a_lo(n), a_hi(n), b_lo(n), b_hi(n);
        for (int i = 0; i < (int)a.size(); i++) a_lo[i] = a[i] % m & 32767, a_hi[i] = a[i] % m >> 15;
        for (int i = 0; i < (int)b.size(); i++) b_lo[i] = b[i] % m & 32767, b_hi[i] = b[i] % m >> 15;
        transform(a_lo, false), transform(a_hi, false), transform(b_lo, false), transform(b_hi, false);

        vector<cd> low(n), mid(n), high(n);
        for (int i = 0; i < n; i++){
            low[i] = a_lo[i] * b_lo[i];
            mid[i] = a_lo[i] * b_hi[i] + a_hi[i] * b_lo[i];
            high[i] = a_hi[i] * b_hi[i];
        }
        transform(low, true), transform(mid, true), transform(high, true);

        vector<long long> res(need);
        for (int i = 0; i < need; i++){
            long long x = llroundl(low[i].real()) % m, y = llroundl(mid[i].real()) % m, z = llroundl(high[i].real()) % m;
            res[i] = (x + (y << 15) % m + (z << 30) % m) % m;
        }
        return res;
    }
}

/// Product of (x + i) for l <= i < r
vector<long long> rising_product(int l, int r, long long m){
    if (r - l == 1) return {l % m, 1 % m};
    int mid = (l + r) / 2;
    return stirling_fft::mod_multiply(rising_product(l, mid, m), rising_product(mid, r, m), m);
}

vector<long long> stirling_first(int n, long long m){
    assert(n >= 0 && 1 <= m && m <= (1LL << 30));
    if (n == 0) return {1 % m};
    return rising_product(0, n, m);
}

vector<long long> stirling_second(int n, long long p){
    assert(n >= 0 && 2 <= p && n < p && p <= (1LL << 30));
    auto pow_mod = [&](long long b, long long e){
        long long res = 1 % p;
        for (b %= p; e; e >>= 1, b = b * b % p){
            if (e & 1) res = res * b % p;
        }
        return res;
    };

    vector<long long> inv_fact(n + 1), x(n + 1), y(n + 1);
    long long fact = 1;
    for (int i = 1; i <= n; i++) fact = fact * i % p;
    inv_fact[n] = pow_mod(fact, p - 2);
    for (int i = n; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % p;

    for (int i = 0; i <= n; i++){
        x[i] = pow_mod(i, n) * inv_fact[i] % p;
        y[i] = i % 2 ? (p - inv_fact[i]) % p : inv_fact[i];
    }
    vector<long long> res = stirling_fft::mod_multiply(x, y, p);
    res.resize(n + 1);
    return res;
}

int main(){
    assert((stirling_first(0, 1000000007) == vector<long long>{1}));
    assert((stirling_first(1, 1000000007) == vector<long long>{0, 1}));
    assert((stirling_first(4, 1000000007) == vector<long long>{0, 6, 11, 6, 1}));
    assert((stirling_first(5, 1000000007) == vector<long long>{0, 24, 50, 35, 10, 1}));
    assert((stirling_first(5, 7) == vector<long long>{0, 3, 1, 0, 3, 1}));
    assert((stirling_first(3, 1) == vector<long long>{0, 0, 0, 0}));

    assert((stirling_second(0, 1000000007) == vector<long long>{1}));
    assert((stirling_second(4, 1000000007) == vector<long long>{0, 1, 7, 6, 1}));
    assert((stirling_second(5, 1000000007) == vector<long long>{0, 1, 15, 25, 10, 1}));
    assert((stirling_second(5, 7) == vector<long long>{0, 1, 1, 4, 3, 1}));
    return 0;
}
