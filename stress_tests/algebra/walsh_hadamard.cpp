#include "../common.h"

#define main library_main
#include "../../code_library/algebra/walsh_hadamard.cpp"
#undef main

/// one shared instance, so every call reuses a scratch buffer left over from a different size
WalshHadamard fwht;

void check_brute(const vector<long long>& a, const vector<long long>& b){
    int n = a.size();
    vector<long long> c_or(n), c_and(n), c_xor(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            c_or[i | j] += a[i] * b[j], c_and[i & j] += a[i] * b[j], c_xor[i ^ j] += a[i] * b[j];
        }
    }

    assert(fwht.or_convolution(a, b) == c_or);
    assert(fwht.and_convolution(a, b) == c_and);
    assert(fwht.xor_convolution(a, b) == c_xor);
}

/// convolving with c * e_k has an O(n) closed form, so it reaches sizes brute force cannot
void check_scaled_unit(int n){
    int k = stress::rand_int(0, n - 1);
    long long c = stress::rand_int(-1000, 1000);
    vector<long long> a(n), b(n);
    for (auto& x : a) x = stress::rand_int(-1000000, 1000000);
    b[k] = c;

    vector<long long> c_or(n), c_and(n), c_xor(n);
    for (int i = 0; i < n; i++) c_or[i | k] += c * a[i], c_and[i & k] += c * a[i], c_xor[i ^ k] += c * a[i];

    assert(fwht.or_convolution(a, b) == c_or);
    assert(fwht.and_convolution(a, b) == c_and);
    assert(fwht.xor_convolution(a, b) == c_xor);
}

long long normalize(long long x, long long mod){
    return (x % mod + mod) % mod;
}

/// inputs may be negative or >= mod, the reference reduces them before multiplying
void check_brute_mod(const vector<long long>& a, const vector<long long>& b, long long mod){
    int n = a.size();
    vector<long long> c_or(n), c_and(n), c_xor(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            long long p = normalize(a[i], mod) * normalize(b[j], mod) % mod;
            c_or[i | j] = (c_or[i | j] + p) % mod, c_and[i & j] = (c_and[i & j] + p) % mod, c_xor[i ^ j] = (c_xor[i ^ j] + p) % mod;
        }
    }

    assert(fwht.or_convolution(a, b, mod) == c_or);
    assert(fwht.and_convolution(a, b, mod) == c_and);
    if (mod % 2) assert(fwht.xor_convolution(a, b, mod) == c_xor);
}

void check_scaled_unit_mod(int n, long long mod){
    int k = stress::rand_int(0, n - 1);
    long long c = stress::rand_int(0, mod - 1);
    vector<long long> a(n), b(n);
    for (auto& x : a) x = stress::rand_int(0, mod - 1);
    b[k] = c;

    vector<long long> c_or(n), c_and(n), c_xor(n);
    for (int i = 0; i < n; i++){
        long long p = c * a[i] % mod;
        c_or[i | k] = (c_or[i | k] + p) % mod, c_and[i & k] = (c_and[i & k] + p) % mod, c_xor[i ^ k] = (c_xor[i ^ k] + p) % mod;
    }

    assert(fwht.or_convolution(a, b, mod) == c_or);
    assert(fwht.and_convolution(a, b, mod) == c_and);
    if (mod % 2) assert(fwht.xor_convolution(a, b, mod) == c_xor);
}

/// all values mod - 1 = -1, so C[k] counts the pairs: 3^popcount(k) for or, 3^(log n - popcount(k)) for and, n for xor
void check_all_max_mod(int log_n, long long mod){
    int n = 1 << log_n;
    vector<long long> a(n, mod - 1), pow3(log_n + 1, 1 % mod);
    for (int i = 1; i <= log_n; i++) pow3[i] = pow3[i - 1] * 3 % mod;

    vector<long long> c_or(n), c_and(n), c_xor(n, n % mod);
    for (int k = 0; k < n; k++) c_or[k] = pow3[__builtin_popcount(k)], c_and[k] = pow3[log_n - __builtin_popcount(k)];

    assert(fwht.or_convolution(a, a, mod) == c_or);
    assert(fwht.and_convolution(a, a, mod) == c_and);
    if (mod % 2) assert(fwht.xor_convolution(a, a, mod) == c_xor);
}

void stress_mod(){
    const vector<long long> mods = {1, 2, 3, 4, 998244353, 1000000007, (1LL << 30), (1LL << 31) - 1};
    for (long long it = 0; it < stress::scaled(600); it++){
        long long mod = it % 3 ? mods[it % mods.size()] : stress::rand_int(1, (1LL << 31) - 1);
        int n = 1 << stress::rand_int(0, it % 20 ? 7 : 10);
        bool wide = stress::rand_int(0, 1);
        vector<long long> a(n), b(n);
        for (auto& x : a) x = wide ? stress::rand_int(-(1LL << 62), 1LL << 62) : stress::rand_int(0, mod - 1);
        for (auto& x : b) x = wide ? stress::rand_int(-(1LL << 62), 1LL << 62) : stress::rand_int(0, mod - 1);
        check_brute_mod(a, b, mod);
    }

    for (long long it = 0; it < stress::scaled(10); it++) check_scaled_unit_mod(1 << stress::rand_int(12, 18), mods[it % mods.size()]);
    check_scaled_unit_mod(1 << 20, 998244353);
    check_scaled_unit_mod(1 << 20, (1LL << 31) - 1);

    check_all_max_mod(20, 998244353);
    check_all_max_mod(20, (1LL << 31) - 1);
    check_all_max_mod(20, 1LL << 30);
    check_all_max_mod(10, 1);
}

/// brute force O(n^2) up to 2^11, closed form for unit vectors up to 2^21, and the n^2 * max|A| * max|B| < 2^62 boundary
/// modular mode: brute force up to 2^10 over small, even and near 2^31 moduli, closed forms at 2^20 with all values mod - 1
int main(){
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = 1 << stress::rand_int(0, it % 20 ? 7 : 11), range = stress::rand_int(0, 1) ? 3 : 100000;
        vector<long long> a(n), b(n);
        for (auto& x : a) x = stress::rand_int(-range, range);
        for (auto& x : b) x = stress::rand_int(-range, range);
        check_brute(a, b);
    }

    for (long long it = 0; it < stress::scaled(20); it++) check_scaled_unit(1 << stress::rand_int(12, 18));
    check_scaled_unit(1 << 21);

    const int n = 1 << 11, limit = (1 << 20) - 1;
    check_brute(vector<long long>(n, limit), vector<long long>(n, -limit));
    check_brute(vector<long long>(n, -limit), vector<long long>(n, -limit));
    for (long long it = 0; it < stress::scaled(4); it++){
        vector<long long> a(n), b(n);
        for (auto& x : a) x = stress::rand_int(0, 1) ? limit : -limit;
        for (auto& x : b) x = stress::rand_int(-limit, limit);
        check_brute(a, b);
    }

    /// near-limit values whose transformed products are close to 2^62 with a wide bit pattern, so lost low bits show
    vector<long long> all_limit(n, limit), near_limit(n, -limit);
    near_limit[0] = -(limit - 1);
    check_brute(all_limit, near_limit);
    for (long long it = 0; it < stress::scaled(4); it++){
        vector<long long> a(n), b(n);
        for (auto& x : a) x = stress::rand_int(limit - 1000, limit);
        for (auto& x : b) x = -stress::rand_int(limit - 1000, limit);
        check_brute(a, b);
    }

    stress_mod();
    return 0;
}
