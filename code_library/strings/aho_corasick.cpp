/***
 *
 * Aho-Corasick Automaton
 * Multi-pattern string matching with failure links
 * Finds all occurrences of multiple pattern strings in a text
 *
 * Complexity:
 *   - insert: O(|s|) plus O(A) to allocate each new trie node, A = MAX_LETTERS = 26
 *   - build: O(A L), L = sum of pattern lengths, A transitions for each of the at most L + 1 nodes
 *   - count: O(text length), counter[] already holds the number of matches ending at each node
 *   - Memory: O(A L) ints for the transition table
 *
 * Usage:
 *   1. Insert all patterns
 *   2. Call build() once, inserting after build() is not supported
 *   3. Query text with count()
 *   4. clear() resets to an empty automaton and keeps vector capacity, so it can be reused for a new pattern set
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED aho_corasick
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
/// END SHARED aho_corasick

int main(){
    auto ac = AhoCorasick();

    ac.insert("hello");
    ac.insert("world");
    ac.insert("lol");
    ac.insert("lol");

    ac.insert("abracadabra");
    ac.insert("abaababbaba");
    ac.insert("aaba");

    ac.build();
    assert(ac.count("helloworldlol") == 4);
    assert(ac.count("abaababbaba") == 2);
    assert(ac.count("aba") == 0);
    assert(ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba") == 7);

    return 0;
}
