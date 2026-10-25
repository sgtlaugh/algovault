#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/josephus_problem.cpp"
#undef main

/// Ids in the order they are killed
vector<int> simulate(int n, long long k){
    vector<int> alive(n), order;
    iota(alive.begin(), alive.end(), 1);
    for (int pos = 0; !alive.empty(); ){
        pos = (pos + (k - 1) % alive.size()) % alive.size();
        order.push_back(alive[pos]);
        alive.erase(alive.begin() + pos);
    }
    return order;
}

/// The O(m) recurrence one circle size at a time, in 128 bits so it is exact for any n and k
long long josephus_linear(long long n, long long k, long long m){
    unsigned __int128 x = n - m;
    for (long long i = n - m + 1; i <= n; i++) x = (x + k) % i;
    return (long long)x + 1;
}

int main(){
    /// Out of range arguments must trip the library's own assert, not some later failure
    auto rejects = [](long long n, long long k, long long m){
        stress::expect_assert("combinatorics/josephus_problem.cpp", [&]{ josephus(n, k, m); });
    };
    rejects(10, 3, 0);
    rejects(10, 3, 11);
    rejects(10, 0, 5);
    rejects(0, 3, 0);

    const long long LIMIT = 4000000000000000000LL;

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 10 ? 30 : 400);
        long long k = it % 3 == 0 ? stress::rand_int(1, 5) : it % 3 == 1 ? stress::rand_int(1, INT_MAX) : stress::rand_int(1, LIMIT);
        auto order = simulate(n, k);
        for (int m = 1; m <= n; m++) assert(josephus(n, k, m) == order[m - 1]);
    }

    /// Too large to simulate: every m up to n against the O(m) recurrence, small k exercises the jumps
    for (long long it = 0; it < stress::scaled(30); it++){
        long long n = stress::rand_int(1, 1000000), m = it % 3 ? n : stress::rand_int(1, n);
        long long k = it % 2 ? stress::rand_int(2, 1000) : stress::rand_int(n, LIMIT);
        assert(josephus(n, k, m) == josephus_linear(n, k, m));
    }

    /// Huge n and k up to the stated limit, the first few thousand kills
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long n = it % 5 == 0 ? LIMIT - stress::rand_int(0, 3) : stress::rand_int(1, LIMIT);
        long long k = it % 4 == 0 ? LIMIT - stress::rand_int(0, 3) : stress::rand_int(1, it % 2 ? 1000 : LIMIT);
        long long m = stress::rand_int(1, min(n, 2000LL));
        assert(josephus(n, k, m) == josephus_linear(n, k, m));
    }

    /// k = 2 survivor has a closed form: n = 2^a + l with l < 2^a gives 2 l + 1, so huge n reach the full jump loop
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long n = it < 10 ? LIMIT - it : stress::rand_int(1, LIMIT);
        long long l = n - (1LL << (63 - __builtin_clzll(n)));
        assert(josephus(n, 2, n) == 2 * l + 1);
    }

    return 0;
}
