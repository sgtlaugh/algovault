#include "../common.h"

#define main library_main
#include "../../code_library/strings/aho_corasick.cpp"
#undef main

string random_word(int len, int alphabet){
    string s(len, 'a');
    for (auto& c : s) c = 'a' + stress::rand_int(0, alphabet - 1);
    return s;
}

long long brute(const vector<string>& patterns, const string& text){
    long long res = 0;
    for (const auto& p : patterns){
        for (size_t i = 0; i + p.size() <= text.size(); i++) res += text.compare(i, p.size(), p) == 0;
    }

    return res;
}

vector<int> brute_each(const vector<string>& patterns, const string& text){
    vector<int> res;
    for (const auto& p : patterns){
        int c = 0;
        for (size_t i = 0; i + p.size() <= text.size(); i++) c += text.compare(i, p.size(), p) == 0;
        res.push_back(c);
    }

    return res;
}

/// Checks the node API against the strings the nodes spell, read off the plain trie before build()
void check_nodes(const vector<string>& patterns, int alphabet){
    AhoCorasick ac;
    for (const auto& p : patterns) ac.insert(p);
    vector<string> spelled(ac.nodes());
    map<string, int> node_of = {{"", 0}};
    for (int u = 0; u < ac.nodes(); u++){
        for (int c = 0; c < 26; c++){
            int v = ac.go[u][c];
            if (v) spelled[v] = spelled[u] + char('a' + c), node_of[spelled[v]] = v;
        }
    }
    assert((int)node_of.size() == ac.nodes());
    for (int v = 0; v < ac.nodes(); v++) assert(ac.ends(v) == count(patterns.begin(), patterns.end(), spelled[v]));
    ac.build();

    /// The node of the longest suffix of s that is a node, the root at worst
    auto longest_suffix_node = [&](const string& s, size_t from){
        for (size_t i = from; i < s.size(); i++){
            if (node_of.count(s.substr(i))) return node_of[s.substr(i)];
        }
        return 0;
    };

    vector<int> position(ac.nodes(), -1);
    assert((int)ac.order.size() == ac.nodes() && ac.order[0] == 0);
    for (int i = 0; i < ac.nodes(); i++) position[ac.order[i]] = i;
    for (int i = 0; i < (int)patterns.size(); i++) assert(ac.terminal[i] == node_of[patterns[i]]);

    for (int v = 0; v < ac.nodes(); v++){
        const string& s = spelled[v];
        assert(position[v] != -1);
        assert(ac.link(v) == longest_suffix_node(s, 1));
        assert(v == 0 || position[ac.link(v)] < position[v]);

        long long suffixes = 0;
        for (const auto& p : patterns) suffixes += p.size() <= s.size() && s.compare(s.size() - p.size(), p.size(), p) == 0;
        assert(ac.matches(v) == suffixes);
        for (int c = 0; c < alphabet; c++) assert(ac.next(v, 'a' + c) == longest_suffix_node(s + char('a' + c), 0));
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        vector<string> patterns(stress::rand_int(0, it % 10 ? 6 : 30));
        int alphabet = stress::rand_int(1, it % 5 ? 3 : 26);
        for (auto& p : patterns) p = random_word(stress::rand_int(1, 6), alphabet);
        check_nodes(patterns, alphabet);
    }

    for (long long it = 0; it < stress::scaled(1000); it++){
        int alphabet = stress::rand_int(1, it % 5 ? 3 : 26);
        vector<string> patterns(stress::rand_int(0, it % 10 ? 8 : 60));
        for (auto& p : patterns) p = random_word(stress::rand_int(1, 6), alphabet);

        AhoCorasick ac;
        for (const auto& p : patterns){
            int v = ac.insert(p), cur = 0;
            for (char c : p) cur = ac.go[cur][c - 'a'];  /// still the plain trie before build()
            assert(cur == v);
        }
        ac.build();
        assert(ac.size() == (int)patterns.size());

        for (int q = 0; q < 10; q++){
            string text = random_word(stress::rand_int(0, it % 10 ? 40 : 2000), alphabet);
            for (auto& c : text) if (stress::rand_int(0, 15) == 0) c = "A{ .\xe9\x80"[stress::rand_int(0, 5)];  /// outside the alphabet, including bytes above 0x7f
            assert(ac.count(text) == brute(patterns, text));
            assert(ac.occurrences(text) == brute_each(patterns, text));
        }
    }

    /// A shifted alphabet: digits, with letters and other bytes in the text acting as separators
    for (long long it = 0; it < stress::scaled(300); it++){
        vector<string> patterns(stress::rand_int(1, 10));
        for (auto& p : patterns){
            p = random_word(stress::rand_int(1, 4), 3);
            for (auto& c : p) c = c - 'a' + '0';
        }

        AhoCorasick<10, '0'> ac;
        for (const auto& p : patterns) ac.insert(p);
        ac.build();
        string text = random_word(stress::rand_int(0, 60), 3);
        for (auto& c : text) c = stress::rand_int(0, 7) ? c - 'a' + '0' : "a/:9\xff"[stress::rand_int(0, 4)];
        assert(ac.count(text) == brute(patterns, text));
        assert(ac.occurrences(text) == brute_each(patterns, text));
    }

    return 0;
}
