/***
 *
 * Aho-Corasick Automaton
 * Multi-pattern string matching with failure links
 * Finds all occurrences of multiple pattern strings in a text
 *
 * Complexity:
 *   - insert: O(|s|) plus O(SIGMA) to allocate each new trie node
 *   - build: O(SIGMA L), L = sum of pattern lengths, SIGMA transitions for each of the at most L + 1 nodes
 *   - count, occurrences: O(text length), plus O(L) for occurrences to push hits up the failure links
 *   - Memory: O(SIGMA L) ints for the transition table
 *
 * AhoCorasick<SIGMA, BASE> ac; pattern characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 *   Text characters outside that range end every match. Remap a sparse alphabet such as ACGT to BASE.. first
 * ac.insert(p): adds pattern p and returns its terminal node, a duplicate pattern is counted once per insert
 * ac.build(): call once after the last insert, inserting after build() is not supported
 * ac.count(t): total number of occurrences of all patterns in t, overlaps included
 * ac.occurrences(t)[i]: number of occurrences of the i-th inserted pattern in t, overlaps included
 * ac.next(v, c): automaton transition for string DP, counter[v] = number of inserted patterns that are suffixes of node v
 * ac.size() is the number of inserted patterns, go.size() the number of nodes, node 0 is the root
 * ac.clear(): resets to an empty automaton and keeps vector capacity, so it can be reused for a new pattern set
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
/// END SHARED aho_corasick

string read_line(){
    static char buf[1 << 17];
    if (!fgets(buf, sizeof(buf), stdin)) return "";
    string s = buf;
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

/// A batch of patterns at most this long keeps the SIGMA-wide transition table near 25 MB
const int BATCH_CHARS = 100000;

int main(){
    int k = stoi(read_line());
    AhoCorasick<58, 'A'> ac;  /// 'A'..'z', the six symbols between 'Z' and 'a' never occur
    while (k--){
        string text = read_line();
        int q = stoi(read_line());
        vector<string> patterns(q);
        for (auto& p : patterns) p = read_line();

        vector<bool> found(q);
        for (int lo = 0, hi; lo < q; lo = hi){
            ac.clear();
            int chars = 0;
            for (hi = lo; hi < q && (hi == lo || chars + (int)patterns[hi].size() <= BATCH_CHARS); hi++){
                chars += patterns[hi].size();
                ac.insert(patterns[hi]);
            }
            ac.build();

            vector<int> hits = ac.occurrences(text);
            /// An empty pattern ends at the root, which counts no hits, yet it is a substring of every text
            for (int i = lo; i < hi; i++) found[i] = hits[i - lo] > 0 || patterns[i].empty();
        }

        for (int i = 0; i < q; i++) puts(found[i] ? "y" : "n");
    }
    return 0;
}
