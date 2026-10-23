/***
 *
 * Aho-Corasick Automaton
 * Counts the occurrences of many patterns in a text in one pass, and walks the automaton for string DP
 *
 * Complexity:
 *   - insert: O(|p|), plus O(SIGMA) for each new node
 *   - build: O(SIGMA L), L = sum of pattern lengths, as there are at most L + 1 nodes
 *   - count: O(|t|)
 *   - occurrences: O(|t| + L)
 *   - memory: O(SIGMA) ints per node
 *
 * AhoCorasick<SIGMA = 26, BASE = 'a'> ac: an empty automaton, pattern characters must lie in [BASE, BASE + SIGMA)
 * ac.build(): finishes the automaton, called once after the last insert
 * ac.clear(): resets to an empty automaton for a new pattern set, keeping vector capacity
 * ac.insert(p): adds pattern p before build() and returns its node, a duplicate pattern counts once per insert
 * ac.count(t): occurrences of all patterns in t, overlaps included
 * ac.ends(v): number of inserted patterns equal to node v
 * ac.link(v): the node of the longest proper suffix of node v, link(0) = 0, after build()
 * ac.matches(v): number of inserted patterns that are suffixes of node v, after build()
 * ac.next(v, c): the node of the longest suffix of node v + c, after build()
 * ac.nodes(): number of nodes, node 0 is the root and each node spells a prefix of a pattern
 * ac.occurrences(t): a vector whose i-th value is the occurrences of the i-th inserted pattern in t
 * ac.size(): number of inserted patterns
 * ac.terminal[i]: the node of the i-th inserted pattern
 * ac.order: every node in BFS order, so link(v) comes before v, after build()
 *
 * A text character outside [BASE, BASE + SIGMA) ends every match, next(v, c) returns the root for it
 * Remap a sparse alphabet such as ACGT to BASE, BASE + 1, ... first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED aho_corasick
template <int SIGMA = 26, char BASE = 'a'>
struct AhoCorasick{
    bool built;
    vector<int> leaf, fail, order, terminal;
    vector<long long> counter;
    vector<string> dictionary;

    /// go[node][char] = trie child before build(), full automaton transition after. 0 doubles as "no child" since the root is never a child
    vector<array<int, SIGMA>> go;

    AhoCorasick(){
        clear();
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

    void clear(){
        go.clear(), dictionary.clear(), order.clear(), terminal.clear();
        fail.clear(), leaf.clear(), counter.clear();

        built = false;
        new_node();
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

    inline int ends(int v){
        return leaf[v];
    }

    inline int link(int v){
        assert(built);
        return fail[v];
    }

    /// Before build() counter is still all zeros, so the assert turns a silent wrong count into an abort
    inline long long matches(int v){
        assert(built);
        return counter[v];
    }

    inline int next(int cur, char ch){
        int x = ch - BASE;
        if (x < 0 || x >= SIGMA) return 0;  /// a letter outside the alphabet ends every match
        return go[cur][x];
    }

    inline int nodes(){
        return go.size();
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

    inline int size(){
        return dictionary.size();
    }

    inline int new_node(){
        leaf.push_back(0);
        counter.push_back(0);
        go.push_back({});
        return go.size() - 1;
    }
};
/// END SHARED aho_corasick

int main(){
    /// Counting seven patterns, "lol" inserted twice, so its hits count twice
    AhoCorasick ac;
    assert(ac.insert("hello") == 5);
    ac.insert("world");
    assert(ac.insert("lol") == 13);
    assert(ac.insert("lol") == 13);

    ac.insert("abracadabra");
    ac.insert("abaababbaba");
    ac.insert("aaba");

    ac.build();
    assert(ac.count("helloworldlol") == 4);
    assert(ac.count("abaababbaba") == 2);
    assert(ac.count("aba") == 0);
    assert(ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba") == 7);
    assert((ac.occurrences("helloworldlolol") == vector<int>{1, 1, 2, 2, 0, 0, 0}));
    assert((ac.occurrences("aabaaba") == vector<int>{0, 0, 0, 0, 0, 0, 2}));

    /// A 4-letter alphabet starting at 'A', where X ends every match
    AhoCorasick<4, 'A'> dna;
    dna.insert("AC");
    dna.insert("CA");
    dna.insert("A");
    dna.build();
    assert(dna.count("ACAXA") == 5);
    assert((dna.occurrences("ACACA") == vector<int>{2, 2, 3}));

    /// Walking the automaton of he, she, his, hers by hand
    AhoCorasick walk;
    for (string p : {"he", "she", "his", "hers"}) walk.insert(p);
    walk.build();
    assert(walk.nodes() == 10);

    int v = 0;
    for (char c : string("ushe")) v = walk.next(v, c);
    assert(v == walk.terminal[1]);
    assert(walk.ends(v) == 1);
    /// she and he both end here
    assert(walk.matches(v) == 2);
    assert(walk.link(v) == walk.terminal[0]);

    return 0;
}
