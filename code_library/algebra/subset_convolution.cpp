/***
 *
 * Subset Convolution and SOS DP
 * Sums or any other fold over all subsets or supersets of each bitmask, sum inverses, subset convolution, exp and log
 *
 * Complexity:
 *   - transforms: O(2^n n) time, in place, and O(2^n n) combine calls for a fold
 *   - subset_convolution: O(2^n n^2) time, 2 * 2^n * (n + 1) ints of memory, 176 MB at n = 20
 *   - subset_inverse: O(2^n n^2) time, as much memory as subset_convolution
 *   - subset_exp, subset_log: O(2^n n^2) time, about 100 MB at n = 20
 *
 * subset_mobius(f, mod = 0): f becomes the g with subset_zeta(g, mod) = f
 * subset_zeta(f, mod = 0): replaces f[S] with the sum of f[T] over every T subset of S, the SOS DP
 * subset_zeta(f, combine): replaces f[S] with combine folded over f[T] for every T subset of S
 * superset_mobius(f, mod = 0): f becomes the g with superset_zeta(g, mod) = f
 * superset_zeta(f, mod = 0): replaces f[S] with the sum of f[T] over every T superset of S
 * superset_zeta(f, combine): replaces f[S] with combine folded over f[T] for every T superset of S
 * subset_convolution(a, b, mod): c[S] = sum of a[T] b[S \ T] over T subset of S, |a| = |b|, 1 <= mod < 2^31
 * subset_exp(f, mod): e[S] = sum over the partitions of S into blocks of the product of f[block], f[0] = 0 (mod mod)
 * subset_inverse(a, mod): the c with subset_convolution(a, c, mod) = {1, 0, ..., 0}, a[0] = 1 (mod mod)
 * subset_log(g, mod): the f with f[0] = 0 and subset_exp(f, mod) = g, g[0] = 1 (mod mod)
 *
 * Notes:
 *   - Every array has size 2^n, n >= 0, index S is the set of bits on in S
 *   - Sums take vector<int> or vector<long long>, folds a vector of any copyable type
 *   - mod = 0 is exact and needs the sum of |f[S]| to fit the element type
 *   - mod >= 1 reduces any input, negatives included, gives results in [0, mod) and must fit the element type
 *   - Nothing is divided, so mod need not be prime
 *   - Connected counts: with f[S] the ways to build the whole structure on S, subset_log(f)[S] counts connected ones
 *   - combine(x, y) must be associative and commutative; each f[T] enters a fold once, so it need not be idempotent
 *   - Requires __int128
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
    assert(mod == 0 || mod - 1 <= (long long)numeric_limits<T>::max());  /// a residue must fit T, or it truncates
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

/// The c with subset_convolution(a, c, mod) = h; one division per mask instead of an inverse then a product, 2x faster
template<typename T>
vector<T> sos_divide(const vector<T>& h, const vector<T>& a, long long mod){
    int size = a.size();
    assert(size > 0 && (size & (size - 1)) == 0 && h.size() == a.size() && mod >= 1 && mod <= INT_MAX);
    assert(sos_reduce_mod(a[0] - 1, mod) == 0);

    int n = __builtin_ctz(size), width = n + 1;
    vector<int> fh(size * width), fa(size * width);
    for (int mask = 0; mask < size; mask++){
        fh[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(h[mask], mod);
        fa[mask * width + __builtin_popcount(mask)] = sos_reduce_mod(a[mask], mod);
    }

    sos_transform(fh, width, mod, false, false);
    sos_transform(fa, width, mod, false, false);

    /// Every mask's rank polynomial of a starts with 1, so dividing by it modulo x^(n + 1) needs no modular inverse
    const unsigned long long square = (unsigned long long)mod * mod;
    for (int mask = 0; mask < size; mask++){
        int* q = &fh[mask * width];
        const int* p = &fa[mask * width];
        int bits = __builtin_popcount(mask);
        for (int k = 1; k <= n; k++){
            unsigned long long sum = 0;
            for (int j = 1; j <= min(k, bits); j++){
                sum += (unsigned long long)p[j] * q[k - j];
                sum = min(sum, sum - square);
            }
            q[k] = (q[k] + mod - sum % mod) % mod;
        }
    }

    sos_transform(fh, width, mod, false, true);

    vector<T> c(size);
    for (int mask = 0; mask < size; mask++) c[mask] = fh[mask * width + __builtin_popcount(mask)];
    return c;
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
    const unsigned long long square = (unsigned long long)mod * mod, inv = ~0ULL / mod;
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
            /// A Barrett reduction of sum < 2^62, 12% faster than % at n = 20
            sum -= (unsigned long long)((unsigned __int128)sum * inv >> 64) * mod;
            x[k] = sum >= (unsigned long long)mod ? sum - mod : sum;
        }
    }

    sos_transform(fa, width, mod, false, true);

    vector<T> c(size);
    for (int mask = 0; mask < size; mask++) c[mask] = fa[mask * width + __builtin_popcount(mask)];
    return c;
}

template<typename T = long long>
vector<T> subset_exp(const vector<T>& f, long long mod){
    int size = f.size();
    assert(size > 0 && (size & (size - 1)) == 0 && mod >= 1 && mod <= INT_MAX && sos_reduce_mod(f[0], mod) == 0);

    /// The block holding the top bit is f[T + bit] for some T, and the rest of the set is any partition below bit
    vector<T> res(size, 0);
    res[0] = 1 % mod;
    for (int bit = 1; bit < size; bit <<= 1){
        vector<T> with_bit(f.begin() + bit, f.begin() + 2 * bit), below(res.begin(), res.begin() + bit);
        vector<T> high = subset_convolution(with_bit, below, mod);
        copy(high.begin(), high.end(), res.begin() + bit);
    }

    return res;
}

template<typename T = long long>
vector<T> subset_inverse(const vector<T>& a, long long mod){
    assert(!a.empty());
    vector<T> one(a.size(), 0);
    one[0] = 1;
    return sos_divide(one, a, mod);
}

template<typename T = long long>
vector<T> subset_log(const vector<T>& g, long long mod){
    int size = g.size();
    assert(size > 0 && (size & (size - 1)) == 0 && mod >= 1 && mod <= INT_MAX && sos_reduce_mod(g[0] - 1, mod) == 0);

    /// The sets below bit are already exp of res, and the sets holding bit are f[T + bit] times that exp
    vector<T> res(size, 0);
    for (int bit = 1; bit < size; bit <<= 1){
        vector<T> with_bit(g.begin() + bit, g.begin() + 2 * bit), below(g.begin(), g.begin() + bit);
        vector<T> high = sos_divide(with_bit, below, mod);
        copy(high.begin(), high.end(), res.begin() + bit);
    }

    return res;
}

int main(){
    /// c[3] = a[0] b[3] + a[1] b[2] + a[2] b[1] + a[3] b[0] = 8 + 14 + 18 + 20
    assert((subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 1000) == vector<long long>{5, 16, 22, 60}));

    /// Partitions of {0, 1} with f = {0, 1, 2, 3}: e[3] = f[3] + f[1] f[2] = 5, and subset_log undoes it
    assert((subset_exp({0, 1, 2, 3}, 1000) == vector<long long>{1, 1, 2, 5}));
    assert((subset_log({1, 1, 2, 5}, 1000) == vector<long long>{0, 1, 2, 3}));

    /// c[3] solves a[0] c[3] + a[1] c[2] + a[2] c[1] + a[3] c[0] = 0, so c[3] = 2 + 2 - 5 = -1
    assert((subset_inverse({1, 1, 2, 5}, 1000) == vector<long long>{1, 999, 998, 999}));

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
