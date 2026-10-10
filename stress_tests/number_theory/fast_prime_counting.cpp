#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/fast_prime_counting.cpp"
#undef main

/// Lucy_Hedgehog's O(n^(3/4)) prime count, an algorithm independent of Meissel-Lehmer
long long lucy(long long n){
    long long r = sqrtl(n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;

    vector<long long> lo(r + 1), hi(r + 1);
    for (long long v = 1; v <= r; v++) lo[v] = v - 1, hi[v] = n / v - 1;

    for (long long p = 2; p <= r; p++){
        if (lo[p] == lo[p - 1]) continue;
        long long below = lo[p - 1], p2 = p * p;
        for (long long i = 1; i <= r && n / i >= p2; i++){
            long long d = i * p <= r ? hi[i * p] : lo[n / (i * p)];
            hi[i] -= d - below;
        }
        for (long long v = r; v >= p2; v--) lo[v] -= lo[v / p] - below;
    }

    return hi[1];
}

/// Tables far smaller than the queries, up to the maxv^2 - 1 limit, and both ends of the MAXP range
template<int MAXN, int MAXP>
void check_small_tables(int maxv){
    PrimeCounter<MAXN, MAXP> pc(maxv);
    long long top = (long long)maxv * maxv - 1;
    assert(pc.lehmer(0) == 0);
    for (long long n : {1LL, 2LL, (long long)maxv - 1, (long long)maxv, top - 1, top}) assert((long long)pc.lehmer(n) == lucy(n));
    for (long long it = 0; it < stress::scaled(30); it++){
        long long n = stress::rand_int(1, top);
        assert((long long)pc.lehmer(n) == lucy(n));
    }
}

int main(){
    PrimeCounter<> pc;
    const long long maxv = pc.maxv;

    /// Around the table boundary and around prime squares and cubes, where phi switches branches
    vector<long long> edges = {maxv - 1, maxv, maxv + 1, 4 * maxv};
    for (long long p : {409LL, 1291LL, 2153LL, 4481LL, 31607LL, 46337LL}){
        for (long long d = -1; d <= 1; d++) edges.push_back(p * p + d), edges.push_back(p * p * p + d);
    }
    for (long long n : edges) if (n <= 5000000000LL) assert((long long)pc.lehmer(n) == lucy(n));

    for (long long it = 0; it < stress::scaled(12); it++){
        long long n = stress::rand_int(it % 3 ? maxv : 1, it % 3 ? 2000000000LL : maxv);
        assert((long long)pc.lehmer(n) == lucy(n));
    }

    for (long long it = 0; it < stress::scaled(1); it++){
        long long n = stress::rand_int(10000000000LL, 30000000000LL);
        assert((long long)pc.lehmer(n) == lucy(n));
    }

    check_small_tables<1, 1>(100);
    check_small_tables<10, 3>(1000);
    check_small_tables<8, 8>(5000);

    return 0;
}
