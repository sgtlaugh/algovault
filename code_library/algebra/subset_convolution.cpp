/***
 *
 * Subset Convolution and SOS DP
 * Zeta and Mobius transforms over subsets and supersets, and subset convolution, modulo mod
 *
 * Complexity: transforms O(2^n n), subset convolution O(2^n n^2) time and O(2^n n) memory
 *
 * Every array has size 2^n, index S is the set of bits on in S
 * Values may be any int, negatives included, and are reduced modulo mod, 1 <= mod < 2^31
 * Nothing is divided, so mod does not have to be prime
 *
 *   subset_zeta(f, mod)              f[S] = sum of f[T] over every T subset of S (SOS DP)
 *   subset_mobius(f, mod)            inverse of subset_zeta
 *   superset_zeta(f, mod)            f[S] = sum of f[T] over every T superset of S
 *   superset_mobius(f, mod)          inverse of superset_zeta
 *   subset_convolution(a, b, mod)    c[S] = sum of a[T] * b[S \ T] over every T subset of S
 *
 * The transforms work in place, subset_convolution allocates 2 * 2^n * (n + 1) ints, 176 MB at n = 20
 *
 * Example:
 *   subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 1000) == {5, 16, 22, 60}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

int sos_add_mod(int a, int b, int mod){
    unsigned sum = (unsigned)a + b;
    return min(sum, sum - mod);
}

int sos_sub_mod(int a, int b, int mod){
    unsigned diff = (unsigned)a - b;
    return min(diff, diff + mod);
}

int sos_reduce_mod(long long x, int mod){
    x %= mod;
    return x < 0 ? x + mod : x;
}

/// Entry (S, k) lives at f[S * width + k], so one pass transforms all width ranked arrays of subset_convolution
void sos_transform(vector<int>& f, int width, int mod, bool superset, bool inverse){
    int size = f.size() / width;
    assert(size > 0 && (size & (size - 1)) == 0 && mod >= 1);

    for (int bit = 1; bit < size; bit <<= 1){
        for (int mask = 0; mask < size; mask++){
            if (!(mask & bit)) continue;

            int* to = &f[(superset ? mask ^ bit : mask) * width];
            const int* from = &f[(superset ? mask : mask ^ bit) * width];
            if (inverse){
                for (int k = 0; k < width; k++) to[k] = sos_sub_mod(to[k], from[k], mod);
            }
            else{
                for (int k = 0; k < width; k++) to[k] = sos_add_mod(to[k], from[k], mod);
            }
        }
    }
}

void sos_reduce_and_transform(vector<int>& f, int mod, bool superset, bool inverse){
    for (auto& x : f) x = sos_reduce_mod(x, mod);
    sos_transform(f, 1, mod, superset, inverse);
}

void subset_zeta(vector<int>& f, int mod){
    sos_reduce_and_transform(f, mod, false, false);
}

void subset_mobius(vector<int>& f, int mod){
    sos_reduce_and_transform(f, mod, false, true);
}

void superset_zeta(vector<int>& f, int mod){
    sos_reduce_and_transform(f, mod, true, false);
}

void superset_mobius(vector<int>& f, int mod){
    sos_reduce_and_transform(f, mod, true, true);
}

vector<int> subset_convolution(const vector<int>& a, const vector<int>& b, int mod){
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

    vector<int> c(size);
    for (int mask = 0; mask < size; mask++) c[mask] = fa[mask * width + __builtin_popcount(mask)];
    return c;
}

int main(){
    assert((subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 1000) == vector<int>{5, 16, 22, 60}));
    assert((subset_convolution({1, 2, 3, 4}, {5, 6, 7, 8}, 7) == vector<int>{5, 2, 1, 4}));
    assert((subset_convolution(vector<int>(8, 1), vector<int>(8, 1), 1000) == vector<int>{1, 2, 2, 4, 2, 4, 4, 8}));
    assert((subset_convolution({7}, {9}, 5) == vector<int>{3}));
    assert((subset_convolution({-1, 4}, {3, -2}, 10) == vector<int>{7, 4}));
    assert((subset_convolution({INT_MIN, INT_MAX}, {1, 1}, 1) == vector<int>{0, 0}));
    assert((subset_convolution({2147483646, 2147483646}, {2147483646, 2147483646}, 2147483647) == vector<int>{1, 2}));
    assert((subset_convolution({2147483646, 2147483646}, {1, 1}, 2147483647) == vector<int>{2147483646, 2147483645}));

    vector<int> f = {1, 2, 3, 4};
    subset_zeta(f, 1000);
    assert((f == vector<int>{1, 3, 4, 10}));
    subset_mobius(f, 1000);
    assert((f == vector<int>{1, 2, 3, 4}));

    superset_zeta(f, 1000);
    assert((f == vector<int>{10, 6, 7, 4}));
    superset_mobius(f, 1000);
    assert((f == vector<int>{1, 2, 3, 4}));

    vector<int> g(8, 1);
    subset_zeta(g, 1000);
    assert((g == vector<int>{1, 2, 2, 4, 2, 4, 4, 8}));
    superset_zeta(g, 5);
    assert((g == vector<int>{2, 3, 3, 2, 3, 2, 2, 3}));

    vector<int> h = {-3, 12, 5, 5};
    subset_zeta(h, 7);
    assert((h == vector<int>{4, 2, 2, 5}));

    return 0;
}
