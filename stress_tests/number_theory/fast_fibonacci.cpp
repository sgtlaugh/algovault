#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/fast_fibonacci.cpp"
#undef main

/// F(n) mod m by squaring the matrix {{1, 1}, {1, 0}}, an independent method
long long matrix_fib(long long n, long long m){
    typedef array<array<__int128, 2>, 2> M;
    auto mul = [&](const M& x, const M& y){
        M z{};
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++) z[i][j] = (x[i][0] * y[0][j] + x[i][1] * y[1][j]) % m;
        }
        return z;
    };
    M res{}, base{};
    res[0][0] = res[1][1] = 1 % m;
    base[0][0] = base[0][1] = base[1][0] = 1 % m;
    for (; n; n >>= 1, base = mul(base, base)){
        if (n & 1) res = mul(res, base);
    }
    return (long long)res[0][1];
}

int main(){
    for (long long m : {1LL, 2LL, 10LL, 1000000007LL, LLONG_MAX, 4611686018427387847LL}){
        long long a = 0, b = 1 % m;
        for (long long n = 0; n <= 20000; n++){
            assert(fibonacci(n, m) == a);
            tie(a, b) = make_pair(b, (long long)(((__int128)a + b) % m));
        }
    }
    long long a = 0, b = 1;
    for (int n = 0; n <= 92; n++){
        assert(fibonacci(n) == a);
        if (n < 92){
            long long next = n + 1 < 92 ? a + b : 0;  /// F(93) would overflow and is never needed
            a = b, b = next;
        }
    }

    for (long long it = 0; it < stress::scaled(50000); it++){
        long long n = it % 3 ? stress::rand_int(0, LLONG_MAX) : LLONG_MAX - stress::rand_int(0, 100);
        long long m = it % 4 == 0 ? LLONG_MAX - stress::rand_int(0, 100) : stress::rand_int(1, LLONG_MAX);
        auto [f, g] = fibonacci_pair(n, m);
        assert(f == matrix_fib(n, m));
        if (n < LLONG_MAX) assert(g == matrix_fib(n + 1, m));
    }
    return 0;
}
