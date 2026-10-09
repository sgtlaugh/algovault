#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sliding_window_aggregation.cpp"
#undef main

const long long MOD = 998244353;

using Linear = pair<long long, long long>;

Linear compose(const Linear& f, const Linear& g){
    return Linear(f.first * g.first % MOD, (f.second * g.first + g.second) % MOD);
}

string concat(const string& a, const string& b){
    return a + b;
}

/// Only a left identity: op(x, "") = "" != x, so the fold must never combine the identity on the right
string newest(const string&, const string& b){
    return b;
}

template <typename T, typename Op>
T naive_fold(const deque<T>& q, T identity, Op op){
    T res = identity;
    for (const T& x : q) res = op(res, x);
    return res;
}

/// Every push/pop pattern up to length 14 with distinct letters, so any misordered or lost element changes the string
/// The newest fold checks that a left identity alone is enough
void exhaustive_patterns(){
    for (int len = 0; len <= 14; len++){
        for (int mask = 0; mask < (1 << len); mask++){
            SlidingWindowAggregation window(string(), concat), last(string(), newest);
            deque<string> ref;
            char next = 'a';
            for (int i = 0; i < len; i++){
                if ((mask >> i & 1) && !ref.empty()) window.pop(), last.pop(), ref.pop_front();
                else{
                    string s(1, next++);
                    window.push(s), last.push(s), ref.push_back(s);
                }
                assert(window.fold() == naive_fold(ref, string(), concat));
                assert(last.fold() == naive_fold(ref, string(), newest));
                assert(window.size() == (int)ref.size() && window.empty() == ref.empty());
            }
        }
    }
}

/// Random interleavings recomputing the fold of non-commutative linear composition mod p, and min/max by scanning
void random_operations(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int ops = stress::rand_int(1, 300), push_bias = stress::rand_int(1, 9);
        long long range = it % 2 ? 3 : 1000000000;
        SlidingWindowAggregation maps(Linear(1, 0), compose);
        MonotonicQueue<long long> low;
        MonotonicQueue<long long, greater<long long>> high;
        deque<Linear> ref_maps;
        deque<long long> ref_values;

        for (int op = 0; op < ops; op++){
            if (!ref_maps.empty() && stress::rand_int(0, 9) >= push_bias){
                maps.pop(), low.pop(), high.pop();
                ref_maps.pop_front(), ref_values.pop_front();
            }
            else{
                Linear f(stress::rand_int(0, MOD - 1), stress::rand_int(0, MOD - 1));
                long long x = stress::rand_int(-range, range);
                maps.push(f), low.push(x), high.push(x);
                ref_maps.push_back(f), ref_values.push_back(x);
            }

            assert(maps.fold() == naive_fold(ref_maps, Linear(1, 0), compose));
            assert(maps.size() == (int)ref_maps.size() && low.size() == (int)ref_values.size());
            assert(low.empty() == ref_values.empty() && high.empty() == ref_values.empty());
            if (ref_values.empty()) continue;
            assert(low.best() == *min_element(ref_values.begin(), ref_values.end()));
            assert(high.best() == *max_element(ref_values.begin(), ref_values.end()));
        }
    }
}

/// Long runs: the min fold must match the monotonic queue, and a periodic scan checks both, which also times the amortized bound
void long_runs(){
    auto take_min = [](long long a, long long b){ return min(a, b); };
    SlidingWindowAggregation window(LLONG_MAX, take_min);
    MonotonicQueue<long long> low;
    deque<long long> ref;

    for (long long op = 0; op < stress::scaled(2000000); op++){
        if (!ref.empty() && stress::rand_int(0, 1)) window.pop(), low.pop(), ref.pop_front();
        else{
            long long x = stress::rand_int(LLONG_MIN + 1, LLONG_MAX - 1);
            window.push(x), low.push(x), ref.push_back(x);
        }

        assert(window.fold() == (ref.empty() ? LLONG_MAX : low.best()));
        if (op % 4096 == 0 && !ref.empty()) assert(low.best() == *min_element(ref.begin(), ref.end()));
    }
}

int main(){
    exhaustive_patterns();
    random_operations();
    long_runs();

    return 0;
}
