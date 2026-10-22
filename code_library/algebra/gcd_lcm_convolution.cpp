/***
 *
 * GCD and LCM Convolution
 * c[k] = sum of a[i] * b[j] over gcd(i, j) = k, or over lcm(i, j) = k, exact or modulo mod
 *
 * Complexity:
 *   - DivisorLattice(n): O(n log log n)
 *   - each convolution or transform: O(m log log m) for vectors of size m + 1
 *
 * DivisorLattice lattice(n): serves any number of calls on vectors of size at most n + 1
 * lattice.gcd_convolution(a, b, mod = 0): c of size max(a.size(), b.size()) with c[k] = sum of a[i] * b[j] over gcd(i, j) = k
 * lattice.lcm_convolution(a, b, mod = 0): the same over lcm(i, j) = k, dropping pairs whose lcm is out of range
 * lattice.multiple_zeta(a, mod = 0): in place, a[d] becomes the sum of a[m] over the multiples m of d
 * lattice.multiple_mobius(a, mod = 0): in place, undoes multiple_zeta
 * lattice.divisor_zeta(a, mod = 0): in place, a[m] becomes the sum of a[d] over the divisors d of m
 * lattice.divisor_mobius(a, mod = 0): in place, undoes divisor_zeta
 *
 * Vectors are indexed by value from 1: the transforms ignore a[0] and the convolutions return c[0] = 0
 * mod = 0 is exact over long long, needing (sum |a[i]|) * (sum |b[j]|) < 2^63 for a convolution and sum |a[i]| < 2^63 for a transform
 * 1 <= mod < 2^62 gives results in [0, mod) for any long long input, negatives included, and mod need not be prime
 * Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DivisorLattice{
    int n;
    vector<int> primes;

    DivisorLattice(int n) : n(n){
        vector<bool> composite(n + 1, false);

        for (int i = 2; i <= n; i++){
            if (composite[i]) continue;
            primes.push_back(i);
            for (long long j = (long long)i * i; j <= n; j += i) composite[j] = true;
        }
    }

    vector<long long> gcd_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, true);
    }

    vector<long long> lcm_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, false);
    }

    void multiple_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) add(a[i], a[i * p], mod);
        }
    }

    void multiple_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) subtract(a[i], a[i * p], mod);
        }
    }

    void divisor_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) add(a[i * p], a[i], mod);
        }
    }

    void divisor_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) subtract(a[i * p], a[i], mod);
        }
    }

    vector<long long> convolve(vector<long long>& a, vector<long long>& b, long long mod, bool by_gcd) const{
        size_t size = max(a.size(), b.size());
        if (!size) return {};
        a.resize(size, 0), b.resize(size, 0);

        if (by_gcd) multiple_zeta(a, mod), multiple_zeta(b, mod);
        else divisor_zeta(a, mod), divisor_zeta(b, mod);

        for (size_t i = 1; i < size; i++){
            a[i] = mod ? (long long)((__int128)a[i] * b[i] % mod) : a[i] * b[i];
        }

        if (by_gcd) multiple_mobius(a, mod);
        else divisor_mobius(a, mod);
        a[0] = 0;
        return a;
    }

    int prepare(vector<long long>& a, long long mod) const{
        assert((long long)a.size() <= (long long)n + 1);
        if (mod) for (auto& x : a) x = (x % mod + mod) % mod;
        return (int)a.size() - 1;
    }

    static void add(long long& x, long long y, long long mod){
        x += y;
        if (mod && x >= mod) x -= mod;
    }

    static void subtract(long long& x, long long y, long long mod){
        x -= y;
        if (mod && x < 0) x += mod;
    }
};

int main(){
    /// gcd and lcm convolutions of a[i] = i and b[j] = j + 4 over 1..4
    const DivisorLattice lattice(4);
    const vector<long long> a = {0, 1, 2, 3, 4}, b = {0, 5, 6, 7, 8};
    assert((lattice.gcd_convolution(a, b) == vector<long long>{0, 155, 52, 21, 32}));  /// c[2] = 2 * 6 + 2 * 8 + 4 * 6
    assert((lattice.lcm_convolution(a, b) == vector<long long>{0, 5, 28, 43, 100}));   /// c[2] = 1 * 6 + 2 * 5 + 2 * 6
    assert((lattice.gcd_convolution(a, b, 7) == vector<long long>{0, 1, 3, 0, 4}));

    /// divisor_zeta of all ones counts divisors, and divisor_mobius undoes it
    vector<long long> f = {0, 1, 1, 1, 1};
    lattice.divisor_zeta(f);
    assert((f == vector<long long>{0, 1, 2, 2, 3}));
    lattice.divisor_mobius(f);
    assert((f == vector<long long>{0, 1, 1, 1, 1}));

    return 0;
}
