#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/fast_prime_sums.cpp"
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

/// Tables far smaller than the queries, up to the maxv^2 - 1 limit
void check_small_tables(int maxv, int maxn, int maxm){
    PrimeSum ps(maxv, maxn, maxm);
    long long top = (long long)maxv * maxv - 1;
    assert(ps.prime_sum(0) == 0);
    for (long long n : {1LL, 2LL, (long long)maxv - 1, (long long)maxv, top - 1, top}) assert(ps.prime_sum(n) == lucy(n));
    for (long long it = 0; it < stress::scaled(30); it++){
        long long n = stress::rand_int(1, top);
        assert(ps.prime_sum(n) == lucy(n));
    }
}

int main(){
    check_small_tables(100, 1, 1);
    check_small_tables(1000, 5, 100);
    check_small_tables(5000, 20, 2000);

    {
        PrimeSum ps(5000000, 20, 200000);   /// the 256 MB preset from the header
        for (long long n : {4999999LL, 5000000LL, (long long)UINT_MAX, UINT_MAX + 1LL}) assert(ps.prime_sum(n) == lucy(n));
        for (long long it = 0; it < stress::scaled(3); it++){
            long long n = stress::rand_int(1, 2000000000LL);
            assert(ps.prime_sum(n) == lucy(n));
        }
    }

    PrimeSum ps;
    const long long maxv = ps.maxv;

    /// The table boundary, both sides of the uint64 / __int128 switch, and around prime squares where phi changes branch
    vector<long long> edges = {maxv - 1, maxv, maxv + 1, UINT_MAX - 1LL, UINT_MAX, UINT_MAX + 1LL, 1LL << 33};
    for (long long p : {1409LL, 4481LL, 46337LL, 65521LL}){
        for (long long d = -1; d <= 1; d++) edges.push_back(p * p + d);
    }
    for (long long n : edges) assert(ps.prime_sum(n) == lucy(n));

    for (long long it = 0; it < stress::scaled(10); it++){
        long long n = stress::rand_int(it % 3 ? maxv : 1, it % 3 ? 2000000000LL : maxv);
        assert(ps.prime_sum(n) == lucy(n));
    }

    for (long long it = 0; it < stress::scaled(1); it++){
        long long n = stress::rand_int(10000000000LL, 30000000000LL);  /// the __int128 instantiation
        assert(ps.prime_sum(n) == lucy(n));
    }

    return 0;
}
