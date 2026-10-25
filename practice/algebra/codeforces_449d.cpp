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
    const int bits = 20, mod = 1000000007;
    int n, x;

    scanf("%d", &n);

    vector<int> f(1 << bits, 0);
    for (int i = 0; i < n; i++){
        scanf("%d", &x);
        f[x]++;
    }

    vector<int> pow2(n + 1, 1);
    for (int i = 1; i <= n; i++) pow2[i] = pow2[i - 1] * 2 % mod;

    superset_zeta(f);
    for (auto& v : f) v = (pow2[v] - 1 + mod) % mod;

    superset_mobius(f, mod);
    printf("%d\n", f[0]);
    return 0;
}
