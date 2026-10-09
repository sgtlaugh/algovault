#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/trie.cpp"
#undef main

template <int BITS>
void stress_binary_trie(long long iterations, int ops, unsigned long long max_value, long long max_copies){
    for (long long it = 0; it < iterations; it++){
        BinaryTrie<BITS> bt;
        map<unsigned long long, long long> ms;
        long long total = 0;
        vector<unsigned long long> pool;
        int pool_size = stress::rand_int(1, 8);
        for (int i = 0; i < pool_size; i++){
            int r = stress::rand_int(0, 3);
            pool.push_back(r == 0 ? 0 : r == 1 ? max_value : (unsigned long long)stress::rand_int(0, max_value >> 1) * 2 + stress::rand_int(0, 1));
        }

        auto random_value = [&](){
            return stress::rand_int(0, 2) ? pool[stress::rand_int(0, pool_size - 1)] : (unsigned long long)stress::rand_int(0, max_value >> 1) * 2 + stress::rand_int(0, 1);
        };

        for (int op = 0; op < ops; op++){
            unsigned long long x = random_value();
            if (stress::rand_int(0, 2) || ms.empty()){
                long long c = stress::rand_int(1, max_copies);
                bt.insert(x, c), ms[x] += c, total += c;
            }
            else{
                auto at = next(ms.begin(), stress::rand_int(0, ms.size() - 1));
                long long c = stress::rand_int(1, at->second);
                bt.erase(at->first, c), total -= c;
                if ((at->second -= c) == 0) ms.erase(at);
            }

            unsigned long long q = random_value();
            assert(bt.size() == total);
            assert(bt.count(q) == (ms.count(q) ? ms[q] : 0));
            if (ms.empty()) continue;

            vector<pair<unsigned long long, long long>> xs;
            for (auto& [y, c] : ms) xs.push_back({q ^ y, c});
            sort(xs.begin(), xs.end());
            assert(bt.min_xor(q) == xs.front().first && bt.max_xor(q) == xs.back().first);

            long long k = stress::rand_int(0, total - 1), seen = 0;
            for (auto& [v, c] : xs){
                if (k < seen + c){
                    assert(bt.kth_xor(q, k) == v);
                    break;
                }
                seen += c;
            }
        }
    }
}

/// Trie prefix and exact counts against a list of inserted words, BinaryTrie xor queries against a std::map multiset
int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        int sigma = stress::rand_int(1, 3);
        Trie<3, 'a'> trie;
        vector<string> words;

        auto random_word = [&](){
            string s;
            for (int i = stress::rand_int(0, 6); i; i--) s += char('a' + stress::rand_int(0, sigma - 1));
            return s;
        };

        for (int op = 0; op < 100; op++){
            if (stress::rand_int(0, 1)){
                string s = random_word();
                trie.insert(s), words.push_back(s);
            }

            string q = stress::rand_int(0, 3) ? random_word() : (words.empty() ? "" : words[stress::rand_int(0, words.size() - 1)]);
            if (stress::rand_int(0, 9) == 0) q.insert(stress::rand_int(0, q.size()), 1, "z`d"[stress::rand_int(0, 2)]);  /// outside [BASE, BASE + SIGMA), 'd' is BASE + SIGMA
            int prefix = 0, exact = 0;
            for (auto& w : words) prefix += w.compare(0, q.size(), q) == 0 && w.size() >= q.size(), exact += w == q;
            assert(trie.count_prefix(q) == prefix && trie.count_word(q) == exact);
        }
    }

    Trie<> big;
    for (int i = 0; i < 100000; i++) big.insert(string(10, 'a' + i % 26));
    assert(big.count_prefix("") == 100000 && big.count_word("aaaaaaaaaa") == 3847 && big.count_prefix("zz") == 3846);

    BinaryTrie<1> one;
    for (int mask = 0; mask < 4; mask++){
        for (int x = 0; x < 2; x++){
            if (mask >> x & 1) one.insert(x);
        }
        for (int q = 0; q < 2 && mask; q++){
            int lo = 2, hi = -1;
            for (int x = 0; x < 2; x++){
                if (mask >> x & 1) lo = min(lo, q ^ x), hi = max(hi, q ^ x);
            }
            assert((int)one.min_xor(q) == lo && (int)one.max_xor(q) == hi);
        }
        for (int x = 0; x < 2; x++){
            if (mask >> x & 1) one.erase(x);
        }
        assert(one.size() == 0);
    }

    stress_binary_trie<3>(stress::scaled(3000), 40, 7, 3);
    stress_binary_trie<10>(stress::scaled(1000), 100, 1023, 1000000);
    stress_binary_trie<4>(stress::scaled(1000), 40, 15, 1000000000000LL);  /// per-node counts far beyond 2^31
    stress_binary_trie<30>(stress::scaled(300), 200, (1ULL << 30) - 1, 1);
    stress_binary_trie<63>(stress::scaled(300), 100, (1ULL << 63) - 1, 2);
    stress_binary_trie<64>(stress::scaled(300), 100, ~0ULL, 2);

    BinaryTrie<> large;
    for (int i = 0; i < 200000; i++) large.insert(i * 5000LL);
    for (int i = 0; i < 200000; i += 2) large.erase(i * 5000LL);
    assert(large.size() == 100000 && large.min_xor(5000) == 0 && large.min_xor(0) == 5000);
    assert(large.max_xor(0) == 199999 * 5000ULL && large.kth_xor(0, 1) == 15000);

    return 0;
}
