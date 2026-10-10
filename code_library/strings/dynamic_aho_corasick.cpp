/***
 *
 * Dynamic Aho-Corasick Automaton
 * Multi-pattern string matching with online pattern insertion
 * Uses binary decomposition over static Aho-Corasick automata
 *
 * Complexity, with N patterns inserted and sigma = MAX_LETTERS:
 *   - insert: O(|pattern| sigma log N) amortized, a pattern is rebuilt into at most log N + 1 automata, each build is O(sigma) per node
 *   - count: O(|text| log N), all floor(log2 N) + 1 slots are scanned even when empty, occurrences come from precomputed counters
 *   - memory: O(sigma) ints per trie node, clear() keeps vector capacity so each slot holds on to its peak size
 *
 * No need to call build() - handled automatically
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/// Static Aho-Corasick (building block for dynamic version)
/// BEGIN COPY aho_corasick from code_library/strings/aho_corasick.cpp
struct AhoCorasick{
    static constexpr int MAX_LETTERS = 26;
    int edge[256];

    vector<int> leaf;
    vector<int> fail;
    vector<long long> counter;
    vector<string> dictionary;

    /// go[node][char] = trie child before build(), full automaton transition after. 0 doubles as "no child" since the root is never a child
    vector<array<int, MAX_LETTERS>> go;

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
        go.clear(), dictionary.clear();
        fail.clear(), leaf.clear(), counter.clear();

        new_node();
        // Map lowercase letters to [0, 25]. Change for different alphabet (digits, uppercase, etc)
        memset(edge, -1, sizeof(edge));
        for (int i = 'a'; i <= 'z'; i++) edge[i] = i - 'a';
    }

    AhoCorasick(){
        clear();
    }

    /// Inserting after build() is not supported
    void insert(const char* str){
        int j, x, cur = 0;

        for (j = 0; str[j] != 0; j++){
            x = edge[(unsigned char)str[j]];
            assert(x >= 0);
            if (!go[cur][x]){
                int next_node = new_node();
                go[cur][x] = next_node;
            }
            cur = go[cur][x];
        }

        leaf[cur]++;
        dictionary.push_back(str);
    }

    void insert(const string& str){
        insert(str.c_str());
    }

    /// Build automaton: compute failure links and fill in missing transitions. Call once after all inserts.
    inline void build(){
        fail.assign(go.size(), 0);
        vector<int> Q = {0};

        for (int i = 0; i < (int)Q.size(); i++){
            int u = Q[i];
            if (u) counter[u] = leaf[u] + counter[fail[u]];
            for (int j = 0; j < MAX_LETTERS; j++){
                int v = go[u][j];
                if (!v) go[u][j] = u ? go[fail[u]][j] : 0;
                else{
                    fail[v] = u ? go[fail[u]][j] : 0;
                    Q.push_back(v);
                }
            }
        }
    }

    inline int next(int cur, char ch){
        int x = edge[(unsigned char)ch];
        if (x < 0) return 0;  /// a letter outside the alphabet ends every match
        return go[cur][x];
    }

    /// Total number of occurrences of all words from dictionary in str
    long long count(const char* str){
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
};
/// END COPY aho_corasick

struct DynamicAhoCorasick{
    /// ar[i] holds 2^i patterns or none, a new top slot is added when every slot is full
    vector<AhoCorasick> ar;

    inline void insert(const char* str){
        // Binary decomposition: find first empty slot
        int i, k = 0;
        for (k = 0; k < (int)ar.size() && ar[k].size(); k++) {}
        if (k == (int)ar.size()) ar.emplace_back();

        // Merge all smaller automata into ar[k]
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
        for (auto& a: ar) res += a.count(str);
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

    return 0;
}
