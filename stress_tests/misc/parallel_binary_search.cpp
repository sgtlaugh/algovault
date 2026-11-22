#include "../common.h"

#define main library_main
#include "../../code_library/misc/parallel_binary_search.cpp"
#undef main

/// Hidden answers; apply and check assert the protocol: updates in order after each reset, checks only read
void check_protocol(int q, int t){
    vector<int> answer(q);
    for (int& a : answer) a = stress::rand_int(0, t + 1);

    int applied = -1, rounds = 0;
    vector<int> res = parallel_binary_search(q, t,
        [&](){ applied = 0, rounds++; },
        [&](int k){
            assert(applied == k && k < t);
            applied++;
        },
        [&](int i){
            assert(applied >= 0 && 0 <= i && i < q);
            return applied >= answer[i];
        });

    assert(res == answer);
    int bound = 1;
    while ((1 << (bound - 1)) < t + 2) bound++;
    assert(rounds <= bound);
}

/// Range adds on a small array: query i asks when the sum over its cells first reaches its need, against a full replay
void check_range_adds(int n, int q, int t){
    vector<int> l(t), r(t), v(t);
    for (int k = 0; k < t; k++){
        l[k] = stress::rand_int(0, n - 1), r[k] = stress::rand_int(l[k], n - 1), v[k] = stress::rand_int(0, 5);
    }
    vector<vector<int>> cells(q);
    vector<long long> need(q);
    for (int i = 0; i < q; i++){
        for (int c = 0; c < n; c++) if (stress::rand_int(0, 2) == 0) cells[i].push_back(c);
        need[i] = stress::rand_int(0, 6 * t + 2);
    }

    vector<long long> a(n);
    vector<int> expected(q, t + 1);
    for (int k = 0; k <= t; k++){
        if (k > 0) for (int c = l[k - 1]; c <= r[k - 1]; c++) a[c] += v[k - 1];
        for (int i = 0; i < q; i++){
            long long sum = 0;
            for (int c : cells[i]) sum += a[c];
            if (expected[i] == t + 1 && sum >= need[i]) expected[i] = k;
        }
    }

    vector<long long> b(n);
    vector<int> res = parallel_binary_search(q, t,
        [&](){ fill(b.begin(), b.end(), 0); },
        [&](int k){ for (int c = l[k]; c <= r[k]; c++) b[c] += v[k]; },
        [&](int i){
            long long sum = 0;
            for (int c : cells[i]) sum += b[c];
            return sum >= need[i];
        });
    assert(res == expected);
}

int main(){
    assert(parallel_binary_search(0, 5, []() {}, [](int) {}, [](int){ return true; }).empty());
    assert((parallel_binary_search(3, 0, []() {}, [](int) {}, [](int i){ return i != 1; }) == vector<int>{0, 1, 0}));

    for (int t = 0; t <= 12; t++){
        for (int q = 0; q <= 6; q++) check_protocol(q, t);
    }
    for (long long it = 0; it < stress::scaled(300); it++){
        check_protocol(stress::rand_int(0, 200), stress::rand_int(0, 1000));
    }
    check_protocol(5000, 100000);

    for (long long it = 0; it < stress::scaled(300); it++){
        check_range_adds(stress::rand_int(1, 8), stress::rand_int(0, 12), stress::rand_int(0, 30));
    }

    stress::expect_assert("misc/parallel_binary_search.cpp", [](){
        parallel_binary_search(-1, 3, []() {}, [](int) {}, [](int){ return true; });
    });
    stress::expect_assert("misc/parallel_binary_search.cpp", [](){
        parallel_binary_search(2, -1, []() {}, [](int) {}, [](int){ return true; });
    });

    return 0;
}
