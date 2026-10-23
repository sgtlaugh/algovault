/***
 *
 * GCD and LCM Convolution
 * c[k] = sum of a[i] * b[j] over gcd(i, j) = k, or over lcm(i, j) = k
 *
 * Complexity: O(n log log n) to build, O(s log log s) per call on vectors of size s
 *
 * DivisorLattice lattice(n): reusable for any vectors of size at most n + 1
 * lattice.divisor_mobius(a, mod = 0): a becomes the f with divisor_zeta(f) = a
 * lattice.divisor_zeta(a, mod = 0): each a[k] becomes the sum of a[d] over the divisors d of k
 * lattice.multiple_mobius(a, mod = 0): a becomes the f with multiple_zeta(f) = a
 * lattice.multiple_zeta(a, mod = 0): each a[k] becomes the sum of a[j] over the multiples j < a.size() of k
 * lattice.gcd_convolution(a, b, mod = 0): c of size max(a.size(), b.size()) by gcd
 * lattice.lcm_convolution(a, b, mod = 0): the same by lcm, pairs with lcm(i, j) >= c.size() dropped
 *
 * Notes:
 *   - Index 0 is ignored and c[0] = 0
 *   - mod = 0 is exact while sum |a[i]| < 2^63 for a transform and (sum |a[i]|) (sum |b[j]|) < 2^63 for a convolution
 *   - Otherwise any 1 <= mod < 2^62 works, prime or not
 *   - Requires __int128
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

    vector<long long> gcd_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, true);
    }

    vector<long long> lcm_convolution(vector<long long> a, vector<long long> b, long long mod = 0) const{
        return convolve(a, b, mod, false);
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
};

int main(){
    /// gcd and lcm convolutions of a[i] = i and b[j] = j + 4 over 1..4
    const DivisorLattice lattice(4);
    const vector<long long> a = {0, 1, 2, 3, 4}, b = {0, 5, 6, 7, 8};
    assert((lattice.gcd_convolution(a, b) == vector<long long>{0, 155, 52, 21, 32}));  /// c[2] = 2 * 6 + 2 * 8 + 4 * 6
    assert((lattice.lcm_convolution(a, b) == vector<long long>{0, 5, 28, 43, 100}));   /// c[2] = 1 * 6 + 2 * 5 + 2 * 6
    assert((lattice.gcd_convolution(a, b, 7) == vector<long long>{0, 1, 3, 0, 4}));

    /// Ordered pairs of values from {2, 4, 3} by gcd: count values each k divides, square, then keep gcd exactly k
    vector<long long> cnt = {0, 0, 1, 1, 1};
    lattice.multiple_zeta(cnt);
    assert((cnt == vector<long long>{0, 3, 2, 1, 1}));
    for (auto& x : cnt) x *= x;
    lattice.multiple_mobius(cnt);
    assert((cnt == vector<long long>{0, 4, 3, 1, 1}));

    /// divisor_zeta of all ones counts divisors, and divisor_mobius undoes it
    vector<long long> f = {0, 1, 1, 1, 1};
    lattice.divisor_zeta(f);
    assert((f == vector<long long>{0, 1, 2, 2, 3}));
    lattice.divisor_mobius(f);
    assert((f == vector<long long>{0, 1, 1, 1, 1}));

    return 0;
}
