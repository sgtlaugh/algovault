/***
 *
 * GCD and LCM Convolution
 * c[k] = sum of a[i] * b[j] over gcd(i, j) = k, or over lcm(i, j) = k, for 1 <= k <= n
 *
 * Complexity: O(n log log n) per call, O(n) memory, plus an O(n log log n) sieve once in the constructor
 *
 * DivisorLattice lattice(n) sieves the primes up to n once, then serves any number of calls on
 * vectors of size at most n + 1, so one lattice sized for the largest input serves all of them
 * Vectors are indexed by value, index 0 is ignored by the transforms and is 0 in the convolution results
 * gcd_convolution(a, b, mod), lcm_convolution(a, b, mod): c of size max(a.size(), b.size()), lcm(i, j) out of range is dropped
 * The four transforms work in place on indices 1..a.size() - 1:
 *     multiple_zeta(a, mod):  a[d] = sum of a[m] over multiples m of d,  multiple_mobius is its inverse
 *     divisor_zeta(a, mod):   a[m] = sum of a[d] over divisors d of m,   divisor_mobius is its inverse
 *
 * mod = 0: exact over long long, the convolutions require (sum |a[i]|) * (sum |b[j]|) < 2^63
 *          and the standalone transforms require sum |a[i]| < 2^63
 * 1 <= mod < 2^62: results in [0, mod), any long long input including negatives, mod need not be prime
 * Requires __int128
 *
 * DivisorLattice lattice(4);
 * lattice.gcd_convolution({0, 1, 2, 3, 4}, {0, 5, 6, 7, 8}) = {0, 155, 52, 21, 32}
 * lattice.lcm_convolution({0, 1, 2, 3, 4}, {0, 5, 6, 7, 8}) = {0, 5, 28, 43, 100}
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

    void divisor_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) subtract(a[i * p], a[i], mod);
        }
    }

    void divisor_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) add(a[i * p], a[i], mod);
        }
    }

    vector<long long> gcd_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, true);
    }

    vector<long long> lcm_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, false);
    }

    void multiple_mobius(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = 1; i <= m / p; i++) subtract(a[i], a[i * p], mod);
        }
    }

    void multiple_zeta(vector<long long>& a, long long mod = 0) const{
        int m = prepare(a, mod);

        for (int p : primes){
            if (p > m) break;
            for (int i = m / p; i >= 1; i--) add(a[i], a[i * p], mod);
        }
    }

private:
    static void add(long long& x, long long y, long long mod){
        x += y;
        if (mod && x >= mod) x -= mod;
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

    static void subtract(long long& x, long long y, long long mod){
        x -= y;
        if (mod && x < 0) x += mod;
    }
};

int main(){
    const DivisorLattice lattice(4);
    const vector<long long> a = {0, 1, 2, 3, 4}, b = {0, 5, 6, 7, 8};
    assert((lattice.gcd_convolution(a, b) == vector<long long>{0, 155, 52, 21, 32}));  /// c[4] = a[4] * b[4] only
    assert((lattice.lcm_convolution(a, b) == vector<long long>{0, 5, 28, 43, 100}));   /// c[2] = 1 * 6 + 2 * 5 + 2 * 6
    assert((lattice.gcd_convolution(a, b, 7) == vector<long long>{0, 1, 3, 0, 4}));    /// the same mod 7

    vector<long long> f = {0, 1, 1, 1, 1};
    lattice.divisor_zeta(f);
    assert((f == vector<long long>{0, 1, 2, 2, 3}));  /// number of divisors
    lattice.divisor_mobius(f);
    assert((f == vector<long long>{0, 1, 1, 1, 1}));
    return 0;
}
