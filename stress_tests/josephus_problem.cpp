#include "common.h"

#define main library_main
#include "../code_library/josephus_problem.cpp"
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

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 10 ? 30 : 400);
        long long k = it % 3 == 0 ? stress::rand_int(1, 5) : it % 3 == 1 ? stress::rand_int(1, INT_MAX) : stress::rand_int(1, 1000000000000000000LL);
        auto order = simulate(n, k);
        for (int m = 1; m <= n; m++){
            if (k <= INT_MAX) assert(josephus1(n, k, m) == order[m - 1]);
            assert(josephus2(n, k, m) == order[m - 1]);
        }
    }

    /// Too large to simulate: the two must agree
    for (long long it = 0; it < stress::scaled(30); it++){
        int n = stress::rand_int(1, 1000000), k = stress::rand_int(it % 2 ? 2 : INT_MAX - 1000, it % 2 ? 1000 : INT_MAX), m = stress::rand_int(1, n);
        assert(josephus1(n, k, m) == josephus2(n, k, m));
    }
    return 0;
}
