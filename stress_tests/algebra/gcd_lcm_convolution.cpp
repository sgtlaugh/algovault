#include "../common.h"

#define main library_main
#include "../../code_library/algebra/gcd_lcm_convolution.cpp"
#undef main

using i128 = __int128;

long long normalize(i128 x, long long mod){
    if (!mod) return (long long)x;
    x %= mod;
    return (long long)(x < 0 ? x + mod : x);
}

/// O(n^2) pair loops over the definition, accumulated in __int128
pair<vector<long long>, vector<long long>> brute(vector<long long> a, vector<long long> b, long long mod){
    size_t size = max(a.size(), b.size());
    a.resize(size, 0), b.resize(size, 0);
    vector<i128> by_gcd(size, 0), by_lcm(size, 0);

    for (size_t i = 1; i < size; i++){
        for (size_t j = 1; j < size; j++){
            i128 product = (i128)a[i] * b[j];
            if (mod) product %= mod;
            by_gcd[gcd(i, j)] += product;
            size_t l = i / gcd(i, j) * j;
            if (l < size) by_lcm[l] += product;
        }
    }

    vector<long long> g(size, 0), l(size, 0);
    for (size_t k = 1; k < size; k++){
        g[k] = normalize(by_gcd[k], mod), l[k] = normalize(by_lcm[k], mod);
    }
    return {g, l};
}

int mobius_brute(int x){
    int mu = 1;
    for (int p = 2; p * p <= x; p++){
        if (x % p) continue;
        x /= p;
        if (x % p == 0) return 0;
        mu = -mu;
    }

    return x > 1 ? -mu : mu;
}

/// Each transform against its defining sum, O(n^2)
void check_transforms(const DivisorLattice& lattice, const vector<long long>& a, long long mod){
    int n = (int)a.size() - 1;
    vector<i128> mz(n + 1, 0), mm(n + 1, 0), dz(n + 1, 0), dm(n + 1, 0);

    for (int d = 1; d <= n; d++){
        for (int m = d; m <= n; m += d){
            mz[d] += a[m], dz[m] += a[d];
            mm[d] += (i128)mobius_brute(m / d) * a[m], dm[m] += (i128)mobius_brute(m / d) * a[d];
        }
    }

    vector<long long> x = a, y = a, z = a, w = a;
    lattice.multiple_zeta(x, mod), lattice.multiple_mobius(y, mod), lattice.divisor_zeta(z, mod), lattice.divisor_mobius(w, mod);
    for (int k = 1; k <= n; k++){
        assert(x[k] == normalize(mz[k], mod) && y[k] == normalize(mm[k], mod));
        assert(z[k] == normalize(dz[k], mod) && w[k] == normalize(dm[k], mod));
    }
}

/// O(n log n) harmonic loops, independent of the prime-by-prime transforms, for n too large for brute
vector<long long> harmonic_gcd(const vector<long long>& a, const vector<long long>& b, long long mod){
    int n = (int)a.size() - 1;
    vector<long long> c(n + 1, 0);

    for (int d = n; d >= 1; d--){
        long long sa = 0, sb = 0;
        for (int m = d; m <= n; m += d) sa = (sa + a[m]) % mod, sb = (sb + b[m]) % mod;
        c[d] = sa * sb % mod;
        for (int m = 2 * d; m <= n; m += d) c[d] = (c[d] - c[m] + mod) % mod;
    }
    return c;
}

vector<long long> harmonic_lcm(const vector<long long>& a, const vector<long long>& b, long long mod){
    int n = (int)a.size() - 1;
    vector<long long> sa(n + 1, 0), sb(n + 1, 0), c(n + 1, 0);

    for (int d = 1; d <= n; d++){
        for (int m = d; m <= n; m += d) sa[m] = (sa[m] + a[d]) % mod, sb[m] = (sb[m] + b[d]) % mod;
    }

    for (int m = 1; m <= n; m++){
        c[m] = (c[m] + sa[m] * sb[m]) % mod;
        for (int k = 2 * m; k <= n; k += m) c[k] = (c[k] - c[m] + mod) % mod;
    }
    return c;
}

vector<long long> random_vector(int size, long long range){
    vector<long long> v(size);
    for (auto& x : v) x = stress::rand_int(-range, range);
    return v;
}

/// Alternates a lattice sized exactly to the input with one shared lattice far larger than it
void check_pair(const DivisorLattice& shared, const vector<long long>& a, const vector<long long>& b, long long mod){
    auto [g, l] = brute(a, b, mod);
    const DivisorLattice exact(max((int)max(a.size(), b.size()) - 1, 0));
    const DivisorLattice& lattice = stress::rand_int(0, 1) ? exact : shared;
    assert(lattice.gcd_convolution(a, b, mod) == g);
    assert(lattice.lcm_convolution(a, b, mod) == l);
}

int main(){
    const DivisorLattice shared(2000);
    const long long mods[] = {0, 998244353, 1, 2, 6, 1000000007, (1LL << 62) - 1, (1LL << 62) - 57};
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long mod = mods[it % 8];
        int sa = it < 200 ? it % 20 : stress::rand_int(0, 80), sb = stress::rand_int(0, 1) ? sa : stress::rand_int(0, 80);
        long long range = mod ? (stress::rand_int(0, 1) ? 3 : LLONG_MAX) : (stress::rand_int(0, 1) ? 3 : 1000000);
        check_pair(shared, random_vector(sa, range), random_vector(sb, range), mod);

        if (it % 4 == 0){
            int size = stress::rand_int(0, 120);
            const DivisorLattice exact(max(size - 1, 0));
            check_transforms(stress::rand_int(0, 1) ? exact : shared, random_vector(size, mod ? LLONG_MAX : 1000000), mod);
        }
    }

    /// The exact-mode bound (sum |a|) * (sum |b|) < 2^63, as close to it as random values get
    for (long long it = 0; it < stress::scaled(40); it++){
        int size = stress::rand_int(2, 1500);
        auto a = random_vector(size, 4000000000LL / size), b = random_vector(size, 4000000000LL / size);
        i128 total_a = 0, total_b = 0;
        for (int i = 1; i < size; i++) total_a += a[i] < 0 ? -a[i] : a[i], total_b += b[i] < 0 ? -b[i] : b[i];
        if (total_a * total_b >= ((i128)1 << 63)) continue;
        check_pair(shared, a, b, 0);
    }

    /// The standalone-transform bound sum |a| < 2^63, (size - 1) * (LLONG_MAX / size) sits just under it
    for (long long it = 0; it < stress::scaled(200); it++){
        int size = stress::rand_int(2, 300);
        check_transforms(shared, random_vector(size, LLONG_MAX / size), 0);
    }

    const long long root = 3037000499LL;
    check_pair(shared, {0, 0, 0, -root}, {0, root, 0, 0}, 0);
    check_pair(shared, {0, 0, 0, 0, 0, 0, root}, {0, 0, 0, root, 0, 0, 0}, 0);

    const int n = 1000000;
    const long long mod = 998244353;
    vector<long long> a(n + 1), b(n + 1);
    for (int i = 1; i <= n; i++) a[i] = stress::rand_int(0, mod - 1), b[i] = stress::rand_int(0, mod - 1);
    const DivisorLattice lattice(n);
    assert(lattice.gcd_convolution(a, b, mod) == harmonic_gcd(a, b, mod));
    assert(lattice.lcm_convolution(a, b, mod) == harmonic_lcm(a, b, mod));

    return 0;
}
