#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/fast_factorial.cpp"
#undef main

bool is_prime(long long n){
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; d++){
        if (n % d == 0) return false;
    }
    return true;
}

long long random_prime(long long lo, long long hi){
    long long p;
    do {
        p = stress::rand_int(lo, hi);
    } while (!is_prime(p));
    return p;
}

/// Every n below p against a running product, so every table entry and every in-block offset is read
void check_all(long long p, long long queries){
    FastFactorial fact(p, queries);
    long long cur = 1;
    for (long long n = 0; n < p; n++){
        if (n) cur = cur * n % p;
        assert(fact(n) == cur);
    }
    assert(fact(p) == 0);
    assert(fact(p + stress::rand_int(0, 1000)) == 0);
}

/// Every (q block)! and a random offset in each block against a running product up to p - 1
void check_blocks(long long p, long long queries){
    FastFactorial fact(p, queries);
    long long cur = 1, next_offset = stress::rand_int(1, fact.block - 1);
    for (long long n = 0; n < p; n++){
        if (n) cur = cur * n % p;
        long long r = n % fact.block;
        if (r == 0) assert(fact(n) == cur);
        if (r == next_offset){
            assert(fact(n) == cur);
            next_offset = stress::rand_int(1, fact.block - 1);
        }
    }
    assert(fact(p - 1) == p - 1);
}

/// Wilson and the reflection n! (p - 1 - n)! = (-1)^(n + 1) need no brute force, so they reach any block near 2^31
void check_identities(long long p, long long queries){
    FastFactorial fact(p, queries);
    long long b = fact.block;
    assert(fact(p - 1) == p - 1);
    assert(fact(LLONG_MAX) == 0);

    long long cur = 1;
    for (long long n = 1; n <= 3 * b + 5; n++){
        cur = cur * n % p;
        if (n % b == 0 || n % 997 == 0) assert(fact(n) == cur);  /// the first entries, built by doubling alone
    }

    for (int it = 0; it < 300; it++){
        long long n = stress::rand_int(0, p - 2);
        if (it % 3 == 0) n = min(p - 2, stress::rand_int(1, (p - 1) / b) * b - stress::rand_int(0, 1));
        long long a = fact(n), r = fact(p - 1 - n), sign = n % 2 ? 1 : p - 1;
        assert(a * r % p == sign);
        assert(fact(n + 1) == a * (n + 1) % p);  /// adjacent values across a block boundary when n = q block - 1
    }
}

int main(){
    for (long long p : {0LL, 1LL, 4LL, 561LL, 2147483645LL, 1LL << 31}){
        stress::expect_assert("combinatorics/fast_factorial.cpp", [p]{ FastFactorial fact(p); });
    }
    stress::expect_assert("combinatorics/fast_factorial.cpp", []{ FastFactorial fact(7, 0); });
    stress::expect_assert("combinatorics/fast_factorial.cpp", []{ FastFactorial(7)(-1); });

    /// Every prime below 1500 covers block 1, where no doubling runs, through block 32
    for (long long p = 2; p < 1500; p++){
        if (is_prime(p)) check_all(p, 1);
    }

    for (long long it = 0; it < stress::scaled(4); it++){
        long long p = random_prime(1500, 30000);
        check_all(p, 1);
    }

    /// Many queries shrink block towards MIN_BLOCK, so the table holds many extension chunks
    check_blocks(4194301, 1);
    check_blocks(4194301, 100000);
    for (long long it = 0; it < stress::scaled(1); it++){
        long long p = random_prime(1000000, 4000000);
        check_blocks(p, stress::rand_int(0, 1) ? 1 : stress::rand_int(1, 1000000));
    }

    /// 2^31 - 1 is the largest prime allowed, its products of two values reach 2^62
    check_identities(2147483647, 1);
    check_identities(1000000007, 1000);
    check_identities(998244353, 1);
    check_identities(random_prime(1LL << 30, (1LL << 31) - 1), stress::rand_int(1, 1000));
    return 0;
}
