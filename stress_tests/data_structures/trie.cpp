#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/trie.cpp"
#undef main

/// Prefix and exact counts against a list of inserted words
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

    return 0;
}
