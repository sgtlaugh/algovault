/***
 *
 * Subset Convolution and SOS DP
 * Zeta and Mobius transforms over subsets and supersets, exact or modulo mod, and subset convolution modulo mod
 *
 * Complexity: transforms O(2^n n), subset convolution O(2^n n^2) time and O(2^n n) memory
 *
 * Every array has size 2^n, index S is the set of bits on in S, values are vector<int> or vector<long long>
 * Nothing is divided, so mod does not have to be prime
 *
 *   subset_zeta(f, mod = 0)          f[S] = sum of f[T] over every T subset of S (SOS DP)
 *   subset_mobius(f, mod = 0)        inverse of subset_zeta
 *   superset_zeta(f, mod = 0)        f[S] = sum of f[T] over every T superset of S
 *   superset_mobius(f, mod = 0)      inverse of superset_zeta
 *   subset_convolution(a, b, mod)    c[S] = sum of a[T] * b[S \ T] over every T subset of S, 1 <= mod < 2^31
 *
 * mod = 0: exact, the transforms require sum |f[S]| to fit the element type
 * mod >= 1: any input values, negatives included, are reduced and results lie in [0, mod), mod must fit the element type
 * The transforms work in place, subset_convolution allocates 2 * 2^n * (n + 1) ints, 176 MB at n = 20
 *
 * Example:
 *   vector<long long> f = {1, 2, 3, 4};
 *   subset_zeta(f);                                    // f == {1, 3, 4, 10}
 *   subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 1000) == {5, 16, 22, 60}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long sos_reduce_mod(long long x, long long mod){
    x %= mod;
    return x < 0 ? x + mod : x;
}

/// Entry (S, k) lives at f[S * width + k], so one pass transforms all width ranked arrays of subset_convolution
template<typename T>
void sos_transform(vector<T>& f, int width, long long mod, bool superset, bool inverse){
    using U = make_unsigned_t<T>;  /// operands lie in [0, mod), so a wrapped sum or difference is undone by one min
    const U umod = mod;
    int size = f.size() / width;
    assert(size > 0 && (size & (size - 1)) == 0 && mod >= 0);

    for (int bit = 1; bit < size; bit <<= 1){
        for (int mask = 0; mask < size; mask++){
            if (!(mask & bit)) continue;

            T* to = &f[(superset ? mask ^ bit : mask) * width];
            const T* from = &f[(superset ? mask : mask ^ bit) * width];
            if (!mod){
                if (inverse) for (int k = 0; k < width; k++) to[k] -= from[k];
                else for (int k = 0; k < width; k++) to[k] += from[k];
            }
            else if (inverse){
                for (int k = 0; k < width; k++){
                    U diff = (U)to[k] - (U)from[k];
                    to[k] = min(diff, diff + umod);
                }
            }
            else{
                for (int k = 0; k < width; k++){
                    U sum = (U)to[k] + (U)from[k];
                    to[k] = min(sum, sum - umod);
                }
            }
        }
    }
}

template<typename T>
void sos_reduce_and_transform(vector<T>& f, long long mod, bool superset, bool inverse){
    if (mod) for (auto& x : f) x = sos_reduce_mod(x, mod);
    sos_transform(f, 1, mod, superset, inverse);
}

template<typename T = long long>
void subset_zeta(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, false, false);
}

template<typename T = long long>
void subset_mobius(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, false, true);
}

template<typename T = long long>
void superset_zeta(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, true, false);
}

template<typename T = long long>
void superset_mobius(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, true, true);
}

template<typename T = long long>
vector<T> subset_convolution(const vector<T>& a, const vector<T>& b, int mod){
    int size = a.size();
    assert(size > 0 && (size & (size - 1)) == 0 && b.size() == a.size() && mod >= 1);

    int n = __builtin_ctz(size), width = n + 1;
    vector<int> fa(size * width), fb(size * width);
    for (int mask = 0; mask < size; mask++){
        fa[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(a[mask], mod);
        fb[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(b[mask], mod);
    }

    sos_transform(fa, width, mod, false, false);
    sos_transform(fb, width, mod, false, false);

    /// Rank k only reads ranks <= k, so going down from n lets the product overwrite fa in place
    /// Ranks above popcount(mask) are still zero after the zeta transform, so the terms using them are skipped
    const unsigned long long square = (unsigned long long)mod * mod;
    for (int mask = 0; mask < size; mask++){
        int* x = &fa[mask * width];
        const int* y = &fb[mask * width];
        int bits = __builtin_popcount(mask);
        for (int k = n; k >= 0; k--){
            unsigned long long sum = 0;
            for (int i = max(0, k - bits); i <= min(k, bits); i++){
                sum += (unsigned long long)x[i] * y[k - i];
                sum = min(sum, sum - square);
            }
            x[k] = sum % mod;
        }
    }

    sos_transform(fa, width, mod, false, true);

    vector<T> c(size);
    for (int mask = 0; mask < size; mask++) c[mask] = fa[mask * width + __builtin_popcount(mask)];
    return c;
}

int main(){
    /// c[3] = a[0] b[3] + a[1] b[2] + a[2] b[1] + a[3] b[0] = 8 + 14 + 18 + 20
    assert((subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 1000) == vector<long long>{5, 16, 22, 60}));
    assert((subset_convolution({-1, 4}, {3, -2}, 10) == vector<long long>{7, 4}));  /// -3 and 2 + 12, reduced mod 10

    vector<int> f = {1, 2, 3, 4};
    subset_zeta(f, 1000);
    assert((f == vector<int>{1, 3, 4, 10}));  /// f[3] = f[0] + f[1] + f[2] + f[3]
    subset_mobius(f, 1000);
    assert((f == vector<int>{1, 2, 3, 4}));

    superset_zeta(f, 1000);
    assert((f == vector<int>{10, 6, 7, 4}));  /// f[1] = f[1] + f[3]
    superset_mobius(f, 1000);
    assert((f == vector<int>{1, 2, 3, 4}));

    vector<long long> big(1 << 20, 1000000000);
    subset_zeta(big);
    assert(big.back() == 1000000000LL << 20 && big[5] == 4000000000LL);
    subset_mobius(big);
    assert(big == vector<long long>(1 << 20, 1000000000));
    return 0;
}
