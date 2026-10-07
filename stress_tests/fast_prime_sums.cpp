#include "common.h"

#define main library_main
#include "../code_library/fast_prime_sums.cpp"
#undef main

/// Lucy_Hedgehog's O(n^(3/4)) prime sum, an algorithm independent of Meissel-Lehmer
__int128 lucy(long long n){
    long long r = sqrtl(n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;

    vector<__int128> lo(r + 1), hi(r + 1);
    auto triangle = [](__int128 v){ return v * (v + 1) / 2 - 1; };
    for (long long v = 1; v <= r; v++) lo[v] = triangle(v), hi[v] = triangle(n / v);

    for (long long p = 2; p <= r; p++){
        if (lo[p] == lo[p - 1]) continue;
        __int128 below = lo[p - 1];
        long long p2 = p * p;
        for (long long i = 1; i <= r && n / i >= p2; i++){
            __int128 d = i * p <= r ? hi[i * p] : lo[n / (i * p)];
            hi[i] -= p * (d - below);
        }
        for (long long v = r; v >= p2; v--) lo[v] -= p * (lo[v / p] - below);
    }
    return hi[1];
}

int main(){
    gen();

    /// The table boundary, both sides of the uint64 / __int128 switch, and around prime squares where phi changes branch
    vector<long long> edges = {MAXV - 1, MAXV, MAXV + 1, UINT_MAX - 1LL, UINT_MAX, UINT_MAX + 1LL, 1LL << 33};
    for (long long p : {1409LL, 4481LL, 46337LL, 65521LL}){
        for (long long d = -1; d <= 1; d++) edges.push_back(p * p + d);
    }
    for (long long n : edges) assert(prime_sum(n) == lucy(n));

    for (long long it = 0; it < stress::scaled(10); it++){
        long long n = stress::rand_int(it % 3 ? MAXV : 1, it % 3 ? 2000000000LL : MAXV);
        assert(prime_sum(n) == lucy(n));
    }
    for (long long it = 0; it < stress::scaled(1); it++){
        long long n = stress::rand_int(10000000000LL, 30000000000LL);  /// the __int128 instantiation
        assert(prime_sum(n) == lucy(n));
    }
    return 0;
}
