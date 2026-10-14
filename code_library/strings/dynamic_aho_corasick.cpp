/***
 *
 * Dynamic Aho-Corasick Automaton
 * Multi-pattern string matching with online pattern insertion
 * Uses binary decomposition over static Aho-Corasick automata
 *
 * Complexity, with N patterns inserted:
 *   - insert: O(|pattern| SIGMA log N) amortized, a pattern is rebuilt into at most log N + 1 automata, each build is O(SIGMA) per node
 *   - count: O(|text| log N), one scan per non-empty slot (popcount(N) of them), occurrences come from precomputed counters
 *   - memory: O(SIGMA) ints per trie node, clear() keeps vector capacity so each slot holds on to its peak size
 *
 * DynamicAhoCorasick<SIGMA, BASE> ac; pattern characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * ac.insert(p): adds pattern p, no build() call needed, ac.count(t): total occurrences of all patterns in t
 * Deletion: keep a second instance for removed patterns, answer = added.count(t) - removed.count(t)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Static Aho-Corasick (building block for dynamic version)
/// BEGIN COPY aho_corasick from code_library/strings/aho_corasick.cpp
template <int SIGMA = 26, char BASE = 'a'>
struct AhoCorasick{
    bool built;
    vector<int> leaf, fail, order, terminal;
    vector<long long> counter;
    vector<string> dictionary;

    /// go[node][char] = trie child before build(), full automaton transition after. 0 doubles as "no child" since the root is never a child
    vector<array<int, SIGMA>> go;

    inline int new_node(){
        leaf.push_back(0);
        counter.push_back(0);
        go.push_back({});
        return go.size() - 1;
    }

    inline int size(){
        return dictionary.size();
    }

    void clear(){
        go.clear(), dictionary.clear(), order.clear(), terminal.clear();
        fail.clear(), leaf.clear(), counter.clear();

        built = false;
        new_node();
    }

    AhoCorasick(){
        clear();
    }

    int insert(const char* str){
        assert(!built);
        int cur = 0;

        for (int j = 0; str[j] != 0; j++){
            int x = str[j] - BASE;
            assert(0 <= x && x < SIGMA);
            if (!go[cur][x]){
                int next_node = new_node();
                go[cur][x] = next_node;
            }
            cur = go[cur][x];
        }

        leaf[cur]++;
        dictionary.push_back(str);
        terminal.push_back(cur);
        return cur;
    }

    int insert(const string& str){
        return insert(str.c_str());
    }

    inline void build(){
        assert(!built);
        built = true;
        fail.assign(go.size(), 0);
        order = {0};

        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            if (u) counter[u] = leaf[u] + counter[fail[u]];
            for (int j = 0; j < SIGMA; j++){
                int v = go[u][j];
                if (!v) go[u][j] = u ? go[fail[u]][j] : 0;
                else{
                    fail[v] = u ? go[fail[u]][j] : 0;
                    order.push_back(v);
                }
            }
        }
    }

    inline int next(int cur, char ch){
        int x = ch - BASE;
        if (x < 0 || x >= SIGMA) return 0;  /// a letter outside the alphabet ends every match
        return go[cur][x];
    }

    inline int nodes(){
        return go.size();
    }

    inline int link(int v){
        assert(built);
        return fail[v];
    }

    inline int ends(int v){
        return leaf[v];
    }

    /// Before build() counter is still all zeros, so the assert turns a silent wrong count into an abort
    inline long long matches(int v){
        assert(built);
        return counter[v];
    }

    long long count(const char* str){
        assert(built);
        long long res = 0;
        for (int j = 0, cur = 0; str[j]; j++){
            cur = next(cur, str[j]);
            res += counter[cur];
        }

        return res;
    }

    long long count(const string& str){
        return count(str.c_str());
    }

    /// Hits of a node also end at every node on its failure chain, so they are pushed up in reverse BFS order
    vector<int> occurrences(const string& text){
        assert(built);
        vector<int> hits(go.size(), 0), res;
        for (int j = 0, cur = 0; j < (int)text.size(); j++){
            cur = next(cur, text[j]);
            hits[cur]++;
        }

        for (int i = (int)order.size() - 1; i > 0; i--) hits[fail[order[i]]] += hits[order[i]];
        for (int v : terminal) res.push_back(hits[v]);
        return res;
    }
};
/// END COPY aho_corasick

template <int SIGMA = 26, char BASE = 'a'>
struct DynamicAhoCorasick{
    /// ar[i] holds 2^i patterns or none, a new top slot is added when every slot is full
    vector<AhoCorasick<SIGMA, BASE>> ar;

    inline void insert(const char* str){
        int i, k;
        for (k = 0; k < (int)ar.size() && ar[k].size(); k++) {}
        if (k == (int)ar.size()) ar.emplace_back();

        ar[k].insert(str);
        for (i = 0; i < k; i++){
            for (const auto& s: ar[i].dictionary){
                ar[k].insert(s);
            }
            ar[i].clear();
        }

        ar[k].build();
    }

    inline void insert(const string& str){
        insert(str.c_str());
    }

    long long count(const char* str){
        long long res = 0;
        for (auto& a: ar) if (a.size()) res += a.count(str);
        return res;
    }

    long long count(const string& str){
        return count(str.c_str());
    }
};

int main(){
    auto ac = DynamicAhoCorasick();

    ac.insert("hello");
    ac.insert("world");

    assert(ac.count("lol") == 0);
    ac.insert("lol");
    assert(ac.count("lol") == 1);
    ac.insert("lol");
    assert(ac.count("lol") == 2);

    ac.insert("abracadabra");
    ac.insert("abaababbaba");
    ac.insert("aaba");

    assert(ac.count("helloworldlol") == 4);
    assert(ac.count("abaababbaba") == 2);
    assert(ac.count("aba") == 0);
    assert(ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba") == 7);

    ac.insert("hello");
    ac.insert("world");
    ac.insert("lol");
    ac.insert("lol");

    ac.insert("a");
    ac.insert("baa");

    assert(ac.count("helloworldlol") == 8);
    assert(ac.count("abaababbaba") == 9);
    assert(ac.count("aba") == 2);
    assert(ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba") == 38);

    DynamicAhoCorasick<26, 'A'> added, removed;
    added.insert("AB"), added.insert("B"), added.insert("AB");
    removed.insert("AB");
    assert(added.count("ABAB") - removed.count("ABAB") == 4);

    return 0;
}
