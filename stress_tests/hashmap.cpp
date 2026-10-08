#include "common.h"

#define main library_main
#include "../code_library/hashmap.cpp"
#undef main

/// Random operations against std::unordered_map, small tables keep probe chains long and wrapping around the end
template <typename K>
void run(long long rounds){
    for (long long round = 0; round < rounds; round++){
        int cap = stress::rand_int(1, round % 10 == 0 ? 2000 : 40), mode = stress::rand_int(0, 3);
        HashMap<K, long long> table(cap);
        unordered_map<K, long long> expected;

        auto random_key = [&]() -> K {
            if (mode == 0) return (K)(stress::rand_int(-cap, cap));                         /// few distinct keys, negatives
            if (mode == 1) return (K)((unsigned long long)stress::rand_int(0, 3 * cap) << (8 * sizeof(K) - 13));  /// only high bits differ, 3 * cap < 2^13
            if (mode == 2 && stress::rand_int(0, 2) == 0) return stress::rand_int(0, 1) ? numeric_limits<K>::min() : numeric_limits<K>::max();
            return (K)stress::rng()();
        };

        for (int op = 0, ops = stress::rand_int(20, 20 + 8 * cap); op < ops; op++){
            K key = random_key();
            long long value = stress::rand_int(-500, 500);
            int kind = stress::rand_int(0, 6);
            bool fits = expected.count(key) || (int)expected.size() < cap;

            if (kind == 0 && fits) table.set(key, value), expected[key] = value;
            else if (kind == 1 && fits) table.add(key, value), expected[key] += value;
            else if (kind <= 3) table.erase(key), expected.erase(key);
            else if (kind == 4 && stress::rand_int(0, 39) == 0) table.clear(), expected.clear();
            else{
                auto it = expected.find(key);
                assert(table.contains(key) == (it != expected.end()));
                assert(table.get(key) == (it == expected.end() ? 0 : it->second));
            }
            assert(table.size() == (int)expected.size());
        }
        for (auto& kv : expected) assert(table.contains(kv.first) && table.get(kv.first) == kv.second);
    }
}

int main(){
    run<int>(stress::scaled(1600));
    run<long long>(stress::scaled(1600));
    run<unsigned int>(stress::scaled(400));
    run<short>(stress::scaled(400));
    return 0;
}
