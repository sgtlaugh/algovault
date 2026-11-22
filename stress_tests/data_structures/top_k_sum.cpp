#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/top_k_sum.cpp"
#undef main

template <typename Compare>
long long brute_sum(vector<long long> values, int k){
    sort(values.begin(), values.end(), Compare());
    long long res = 0;
    for (int i = 0; i < min(k, (int)values.size()); i++) res += values[i];

    return res;
}

/// Random insert, erase and set_k with few distinct values, so ties at the k-th position are common
template <typename Compare>
void random_run(int ops, long long lo, long long hi, int max_k){
    int k = stress::rand_int(0, max_k);
    TopKSum<long long, Compare> s(k);
    vector<long long> values;

    for (int it = 0; it < ops; it++){
        int op = stress::rand_int(0, 9);
        if (op < 5 || values.empty()){
            long long x = stress::rand_int(lo, hi);
            s.insert(x);
            values.push_back(x);
        }
        else if (op < 8){
            int i = stress::rand_int(0, (int)values.size() - 1);
            s.erase(values[i]);
            values.erase(values.begin() + i);
        }
        else{
            k = op == 8 ? stress::rand_int(0, max_k) : max(0, k + (int)stress::rand_int(-1, 1));
            s.set_k(k);
        }

        assert(s.size() == (int)values.size());
        assert(s.sum == brute_sum<Compare>(values, k));
    }
}

int main(){
    TopKSum empty(0);
    assert(empty.sum == 0 && empty.size() == 0);

    /// Only orders whose ties are equal values, as the header requires: an |x| order breaks erase by design
    for (int it = 0; it < stress::scaled(3000); it++){
        int max_k = stress::rand_int(0, 8);
        long long range = stress::rand_int(0, 6);
        random_run<greater<long long>>(60, -range, range, max_k);
        random_run<less<long long>>(60, -range, range, max_k);
    }

    for (int it = 0; it < stress::scaled(30); it++){
        random_run<greater<long long>>(2000, -(long long)1e12, (long long)1e12, 300);
        random_run<less<long long>>(2000, -(long long)1e12, (long long)1e12, 300);
    }

    /// Every value erased then reinserted, k jumping past the size and back to 0
    TopKSum<int, less<int>> low(3);
    for (int round = 0; round < 3; round++){
        for (int x = 1; x <= 6; x++) low.insert(x);
        assert(low.sum == 6);
        low.set_k(100);
        assert(low.sum == 21);
        low.set_k(0);
        assert(low.sum == 0);
        low.set_k(3);
        for (int x = 1; x <= 6; x++) low.erase(x);
        assert(low.sum == 0 && low.size() == 0);
    }

    stress::expect_assert("data_structures/top_k_sum.cpp", []{ TopKSum s(2); s.erase(1); });
    stress::expect_assert("data_structures/top_k_sum.cpp", []{ TopKSum s(-1); });
    stress::expect_assert("data_structures/top_k_sum.cpp", []{ TopKSum s(1); s.set_k(-2); });

    return 0;
}
