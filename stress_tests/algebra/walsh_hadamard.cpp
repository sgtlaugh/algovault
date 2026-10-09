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

/// brute force O(n^2) up to 2^11, closed form for unit vectors up to 2^21, and the n^2 * max|A| * max|B| < 2^62 boundary
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

    return 0;
}
