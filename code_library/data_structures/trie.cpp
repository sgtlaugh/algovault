/***
 *
 * Trie
 * Prefix tree over strings, counting how many inserted words pass through or end at each node
 *
 * Complexity: O(|s|) per operation, O(total length * SIGMA) memory
 *
 * Trie<SIGMA, BASE> trie; characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * trie.insert(s): adds one copy of s, duplicates count separately
 * trie.count_prefix(p): inserted words that start with p (every word for the empty prefix)
 * trie.count_word(s): copies of exactly s
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int SIGMA = 26, char BASE = 'a'>
struct Trie{
    vector<array<int, SIGMA>> next;
    vector<int> pass, ends;

    Trie(){
        new_node();
    }

    int new_node(){
        next.push_back({});
        next.back().fill(-1);
        pass.push_back(0), ends.push_back(0);
        return next.size() - 1;
    }

    void insert(const string& s){
        int v = 0;
        pass[0]++;
        for (char ch : s){
            int c = ch - BASE;
            assert(0 <= c && c < SIGMA);
            if (next[v][c] == -1){
                int u = new_node();
                next[v][c] = u;
            }
            v = next[v][c];
            pass[v]++;
        }
        ends[v]++;
    }

    /// Node reached by s, or -1 if no inserted word has s as a prefix
    int find(const string& s) const{
        int v = 0;
        for (char ch : s){
            int c = ch - BASE;
            if (c < 0 || c >= SIGMA || next[v][c] == -1) return -1;
            v = next[v][c];
        }
        return v;
    }

    int count_prefix(const string& p) const{
        int v = find(p);
        return v == -1 ? 0 : pass[v];
    }

    int count_word(const string& s) const{
        int v = find(s);
        return v == -1 ? 0 : ends[v];
    }
};

int main(){
    Trie<> trie;
    for (string s : {"apple", "app", "apply", "banana", "app"}) trie.insert(s);
    assert(trie.count_prefix("app") == 4 && trie.count_prefix("appl") == 2 && trie.count_prefix("") == 5);
    assert(trie.count_prefix("b") == 1 && trie.count_prefix("c") == 0 && trie.count_prefix("applex") == 0);
    assert(trie.count_word("app") == 2 && trie.count_word("appl") == 0 && trie.count_word("banana") == 1);
    assert(trie.count_word("Apple") == 0);

    trie.insert("");
    assert(trie.count_word("") == 1 && trie.count_prefix("") == 6);

    Trie<2, '0'> bits;
    bits.insert("0101"), bits.insert("0110");
    assert(bits.count_prefix("01") == 2 && bits.count_prefix("010") == 1 && bits.count_word("0110") == 1);
    return 0;
}
