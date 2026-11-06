/***
 *
 * Subset Convolution and SOS DP
 * Sums or any other fold over all subsets or supersets of each bitmask, sum inverses, and subset convolution
 *
 * Complexity:
 *   - transforms: O(2^n n) time, in place, and O(2^n n) combine calls for a fold
 *   - subset_convolution: O(2^n n^2) time, 2 * 2^n * (n + 1) ints of memory, 176 MB at n = 20
 *
 * subset_mobius(f, mod = 0): f becomes the g with subset_zeta(g, mod) = f
 * subset_zeta(f, mod = 0): replaces f[S] with the sum of f[T] over every T subset of S, the SOS DP
 * subset_zeta(f, combine): replaces f[S] with combine folded over f[T] for every T subset of S
 * superset_mobius(f, mod = 0): f becomes the g with superset_zeta(g, mod) = f
 * superset_zeta(f, mod = 0): replaces f[S] with the sum of f[T] over every T superset of S
 * superset_zeta(f, combine): replaces f[S] with combine folded over f[T] for every T superset of S
 * subset_convolution(a, b, mod): c[S] = sum of a[T] b[S \ T] over T subset of S, |a| = |b|, 1 <= mod < 2^31
 *
 * Notes:
 *   - Every array has size 2^n, n >= 0, index S is the set of bits on in S
 *   - Sums take vector<int> or vector<long long>, folds a vector of any copyable type
 *   - mod = 0 is exact and needs the sum of |f[S]| to fit the element type
 *   - mod >= 1 reduces any input, negatives included, gives results in [0, mod) and must fit the element type
 *   - Nothing is divided, so mod need not be prime
 *   - combine(x, y) must be associative and commutative; each f[T] enters a fold once, so it need not be idempotent
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

    /// Masks start..start + bit - 1 lack bit and are contiguous, their entries pair with the next bit * width entries
    for (int bit = 1; bit < size; bit <<= 1){
        int len = bit * width;
        for (int start = 0; start < size; start += 2 * bit){
            T* without_bit = &f[start * width];
            T* to = superset ? without_bit : without_bit + len;
            const T* from = superset ? without_bit + len : without_bit;

            if (!mod){
                if (inverse) for (int i = 0; i < len; i++) to[i] -= from[i];
                else for (int i = 0; i < len; i++) to[i] += from[i];
            }
            else if (inverse){
                for (int i = 0; i < len; i++){
                    U diff = (U)to[i] - (U)from[i];
                    to[i] = min(diff, diff + umod);
                }
            }
            else{
                for (int i = 0; i < len; i++){
                    U sum = (U)to[i] + (U)from[i];
                    to[i] = min(sum, sum - umod);
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

template<typename T, typename Combine>
void sos_fold(vector<T>& f, Combine combine, bool superset){
    int size = f.size();
    assert(size > 0 && (size & (size - 1)) == 0);

    for (int bit = 1; bit < size; bit <<= 1){
        for (int start = 0; start < size; start += 2 * bit){
            for (int i = start; i < start + bit; i++){
                if (superset) f[i] = combine(f[i], f[i + bit]);
                else f[i + bit] = combine(f[i + bit], f[i]);
            }
        }
    }
}

template<typename T = long long>
void subset_mobius(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, false, true);
}

template<typename T = long long>
void subset_zeta(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, false, false);
}

/// Without the enable_if, subset_zeta(f, 7) would pick this overload with Combine = int over the long long mod one
template<typename T, typename Combine, enable_if_t<is_invocable_v<Combine&, const T&, const T&>, int> = 0>
void subset_zeta(vector<T>& f, Combine combine){
    sos_fold(f, combine, false);
}

template<typename T = long long>
void superset_mobius(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, true, true);
}

template<typename T = long long>
void superset_zeta(vector<T>& f, long long mod = 0){
    sos_reduce_and_transform(f, mod, true, false);
}

template<typename T, typename Combine, enable_if_t<is_invocable_v<Combine&, const T&, const T&>, int> = 0>
void superset_zeta(vector<T>& f, Combine combine){
    sos_fold(f, combine, true);
}

template<typename T = long long>
vector<T> subset_convolution(const vector<T>& a, const vector<T>& b, long long mod){
    int size = a.size();

    assert(size > 0 && (size & (size - 1)) == 0 && b.size() == a.size() && mod >= 1 && mod <= INT_MAX);

    int n = __builtin_ctz(size), width = n + 1;
    vector<int> fa(size * width), fb(size * width);
    for (int mask = 0; mask < size; mask++){
        fa[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(a[mask], mod);
        fb[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(b[mask], mod);
    }

    sos_transform(fa, width, mod, false, false);
    sos_transform(fb, width, mod, false, false);

    /// Rank k reads ranks <= k only, so going down from n can overwrite fa, and ranks above popcount(mask) are zero
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

    /// Exact sums over subsets and supersets, and back
    vector<long long> f = {1, 2, 3, 4};
    subset_zeta(f);
    assert((f == vector<long long>{1, 3, 4, 10}));  /// f[3] = f[0] + f[1] + f[2] + f[3]
    subset_mobius(f);
    assert((f == vector<long long>{1, 2, 3, 4}));

    superset_zeta(f);
    assert((f == vector<long long>{10, 6, 7, 4}));  /// f[1] = f[1] + f[3]
    superset_mobius(f);
    assert((f == vector<long long>{1, 2, 3, 4}));

    /// Modulo 7, -1 is reduced to 6 first
    vector<int> g = {-1, 2, 3, 4};
    subset_zeta(g, 7);
    assert((g == vector<int>{6, 1, 2, 1}));

    /// Maximum over subsets and supersets with a lambda in place of the sum
    vector<int> h = {5, 2, 7, 1};
    subset_zeta(h, [](int x, int y){ return max(x, y); });
    assert((h == vector<int>{5, 5, 7, 7}));

    h = {5, 2, 7, 1};
    superset_zeta(h, [](int x, int y){ return max(x, y); });
    assert((h == vector<int>{7, 2, 7, 1}));

    return 0;
}
