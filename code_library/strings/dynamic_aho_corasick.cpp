/***
 *
 * Dynamic Aho-Corasick Automaton
 * Multi-pattern string matching with online pattern insertion
 * Uses binary decomposition over static Aho-Corasick automata
 *
 * Insert: O(|pattern| * log N) amortized
 * Query:  O(|text| * log N + occurrences * log N)
 *
 * No need to call build() - handled automatically
 * Holds at most 2^MAX_LOG - 1 patterns
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAX_LOG      20
#define MAX_LETTERS  26

using namespace std;

/// Static Aho-Corasick (building block for dynamic version)
struct AhoCorasick{
    int edge[256];

    vector<int> leaf;
    vector<int> fail;
    vector<long long> counter;
    vector<string> dictionary;

    /// go[node][char] = trie child before build(), full automaton transition after. 0 doubles as "no child" since the root is never a child
    vector<array<int, MAX_LETTERS>> go;

    inline int node(){
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

        node();
        // Map lowercase letters to [0, 25]. Change for different alphabet (digits, uppercase, etc)
        memset(edge, -1, sizeof(edge));
        for (int i = 'a'; i <= 'z'; i++) edge[i] = i - 'a';
    }

    AhoCorasick(){
        clear();
    }

    /// Inserting after build() is not supported, DynamicAhoCorasick rebuilds a fresh automaton instead
    inline void insert(const char* str){
        int j, x, cur = 0;

        for (j = 0; str[j] != 0; j++){
            x = edge[(unsigned char)str[j]];
            assert(x >= 0);
            if (!go[cur][x]){
                int next_node = node();
                go[cur][x] = next_node;
            }
            cur = go[cur][x];
        }

        leaf[cur]++;
        dictionary.push_back(str);
    }

    inline void insert(const string& str){
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

    /// total number of occurrences of all words from dictionary in str
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

struct DynamicAhoCorasick{
    AhoCorasick ar[MAX_LOG];

    inline void insert(const char* str){
        // Binary decomposition: find first empty slot
        int i, k = 0;
        for (k = 0; k < MAX_LOG && ar[k].size(); k++) {}
        assert(k < MAX_LOG);

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
        for (int i = 0; i < MAX_LOG; i++) res += ar[i].count(str);
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
