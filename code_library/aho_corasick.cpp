/***
 *
 * Aho-Corasick Automaton
 * Multi-pattern string matching with failure links
 * Finds all occurrences of multiple pattern strings in a text
 *
 * Build: O(sum of pattern lengths)
 * Query: O(text length + total occurrences)
 *
 * Usage:
 *   1. Insert all patterns
 *   2. Call build() once
 *   3. Query text with count()
 *
***/

#include <bits/stdtr1c++.h>

#define MAX_LETTERS 26

using namespace std;

struct AhoCorasick{
    int id, edge[256];

	vector<int> leaf;
    vector<int> fail;
    vector<long long> counter;
    vector<string> dictionary;

    vector<vector<int>> dp;      // dp[node][char] = precomputed next state (O(1) transitions)
    vector<map<char, int>> trie; // Trie structure: trie[node][char] = child node

    inline int new_node(){
        leaf.push_back(0);
        counter.push_back(0);
        dp.push_back(vector<int>(MAX_LETTERS, 0));
        trie.push_back(map<char, int>());
        return id++;
    }

    inline int size(){
        return dictionary.size();
    }

    AhoCorasick(){
        trie.clear(), dictionary.clear();
        dp.clear(), fail.clear(), leaf.clear(), counter.clear();
        id = 0, new_node();

        // Map lowercase letters to [0, 25]. Change for different alphabet (digits, uppercase, etc)
        for (int i = 'a'; i <= 'z'; i++){
            edge[i] = i - 'a';
        }
    }

    void insert(const char* str){
        int j, x, cur = 0;
        for (j = 0; str[j] != 0; j++){
            x = edge[(int)str[j]];
            if (!trie[cur].count(x)){
                int next_node = new_node();
                trie[cur][x] = next_node;
            }
            cur = trie[cur][x];
        }

        leaf[cur]++;
        dictionary.push_back(str);
    }

    void insert(const string& str){
        insert(str.c_str());
    }

    /// Build automaton: compute failure links and precompute transitions. Call once after all inserts.
    inline void build(){
        vector<pair<int, pair<int, int>>> Q;
        fail.resize(id, 0);
        Q.push_back({0, {0, 0}});

        for (int i = 0; i < id; i++){
            for (int j = 0; j < MAX_LETTERS; j++){
                dp[i][j] = i;
            }
        }

        for(int i = 0; i < (int)Q.size(); i++){
            int u = Q[i].first;
            int p = Q[i].second.first;
            char c = Q[i].second.second;
            for(auto& it: trie[u]){
                Q.push_back({it.second, {u, it.first}});
            }

            if (u){
                int f = fail[p];
                while (f && !trie[f].count(c)) f = fail[f];
                if(!trie[f].count(c) || trie[f][c] == u) fail[u] = 0;
                else fail[u] = trie[f][c];
                counter[u] = leaf[u] + counter[fail[u]];

                for (int j = 0; j < MAX_LETTERS; j++){
                    if (u && !trie[u].count(j)){
                        dp[u][j] = dp[fail[u]][j];
                    }
                }
            }
        }
    }

    inline int next(int cur, char ch){
        int x = edge[(int)ch];
        cur = dp[cur][x];
        if (trie[cur].count(x)) cur = trie[cur][x];
        return cur;
    }

    /// Total number of occurrences of all words from dictionary in str
    long long count(const char* str){
        long long res = 0;
        for (int j = 0, cur = 0; str[j] && id > 1; j++){
            cur = next(cur, str[j]);
            res += counter[cur];
        }

        return res;
    }

    long long count(const string& str){
        return count(str.c_str());
    }
};


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
