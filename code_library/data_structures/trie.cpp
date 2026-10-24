/***
 *
 * Trie and Binary Trie
 * Prefix tree that counts inserted words by prefix, and a binary trie for xor queries over a multiset of integers
 *
 * Complexity:
 *   - Trie: O(|s|) per operation, O(L SIGMA) memory for L inserted characters in total
 *   - BinaryTrie: O(BITS) per operation, O(v BITS) memory for v distinct inserted values
 *
 * Trie<SIGMA = 26, BASE = 'a'> trie: no words yet, characters must lie in [BASE, BASE + SIGMA)
 * Trie<> trie(max_nodes): the same with memory for max_nodes >= 1 nodes reserved up front, more nodes still fit
 * trie.insert(s): adds one copy of s, duplicates count separately
 * trie.child(v, ch): node one character ch below node v, -1 if there is none or ch is outside the alphabet, v != -1
 * trie.count_prefix(p): inserted words that start with p, every word for the empty prefix
 * trie.count_word(s): copies of exactly s
 * trie.find(s): node of the prefix s, -1 if no inserted word starts with s, the root 0 for the empty prefix
 * trie.pass[v]: inserted words that start with node v's prefix
 * trie.ends[v]: copies of exactly node v's prefix
 *
 * BinaryTrie<BITS = 30> bt: an empty multiset of unsigned long long values in [0, 2^BITS), 1 <= BITS <= 64
 * bt.erase(x, c = 1): removes c >= 1 copies of x, count(x) >= c
 * bt.insert(x, c = 1): adds c >= 1 copies of x
 * bt.count(x): copies of x
 * bt.kth_xor(x, k): k-th smallest x ^ y counting copies, 0-based, 0 <= k < size()
 * bt.max_xor(x): max of x ^ y over the y in the multiset, size() > 0
 * bt.min_xor(x): min of x ^ y over the y in the multiset, size() > 0
 * bt.size(): copies of all values
 *
 * Notes:
 *   - A trie holds at most 1 + total inserted length nodes, the natural max_nodes, of 4 * SIGMA + 8 bytes each
 *   - Without max_nodes the vector doubling briefly holds up to twice the memory of the nodes in use
 *   - BinaryTrie keeps erased nodes allocated, so memory grows with insertions of new values, not with live size
 *   - A one-off total over all words' shared prefixes needs no trie: sort the words and sum their adjacent common prefixes
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

    explicit Trie(int max_nodes){
        assert(max_nodes >= 1);
        next.reserve(max_nodes), pass.reserve(max_nodes), ends.reserve(max_nodes);
        new_node();
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

    int child(int v, char ch) const{
        int c = ch - BASE;
        return c < 0 || c >= SIGMA ? -1 : next[v][c];
    }

    /// Deleted so that child(v, s[i] - 'a') fails to compile instead of reading the int as a character
    int child(int v, int c) const = delete;

    int count_prefix(const string& p) const{
        int v = find(p);
        return v == -1 ? 0 : pass[v];
    }

    int count_word(const string& s) const{
        int v = find(s);
        return v == -1 ? 0 : ends[v];
    }

    int find(const string& s) const{
        int v = 0;
        for (int i = 0; i < (int)s.size() && v != -1; i++) v = child(v, s[i]);
        return v;
    }

    int new_node(){
        next.push_back({});
        next.back().fill(-1);
        pass.push_back(0), ends.push_back(0);
        return next.size() - 1;
    }
};

template <int BITS = 30>
struct BinaryTrie{
    static_assert(1 <= BITS && BITS <= 64, "BITS must lie in [1, 64]");

    vector<array<int, 2>> next;
    vector<long long> cnt;

    BinaryTrie(){
        new_node();
    }

    void erase(unsigned long long x, long long c = 1){
        assert(c >= 1 && count(x) >= c);
        int v = 0;
        cnt[0] -= c;
        for (int b = BITS - 1; b >= 0; b--){
            v = next[v][x >> b & 1];
            cnt[v] -= c;
        }
    }

    void insert(unsigned long long x, long long c = 1){
        check_range(x);
        assert(c >= 1);
        int v = 0;
        cnt[0] += c;
        for (int b = BITS - 1; b >= 0; b--){
            int bit = x >> b & 1;
            if (next[v][bit] == -1){
                int u = new_node();
                next[v][bit] = u;
            }
            v = next[v][bit];
            cnt[v] += c;
        }
    }

    long long count(unsigned long long x) const{
        check_range(x);
        int v = 0;
        for (int b = BITS - 1; b >= 0 && v != -1; b--) v = next[v][x >> b & 1];
        return weight(v);
    }

    unsigned long long kth_xor(unsigned long long x, long long k) const{
        check_range(x);
        assert(0 <= k && k < size());
        int v = 0;
        unsigned long long res = 0;

        for (int b = BITS - 1; b >= 0; b--){
            int bit = x >> b & 1;
            int same = next[v][bit];
            if (k < weight(same)){
                v = same;
                continue;
            }
            k -= weight(same);
            v = next[v][bit ^ 1];
            res |= 1ULL << b;
        }

        return res;
    }

    unsigned long long max_xor(unsigned long long x) const{
        return kth_xor(x, size() - 1);
    }

    unsigned long long min_xor(unsigned long long x) const{
        return kth_xor(x, 0);
    }

    long long size() const{
        return cnt[0];
    }

    int new_node(){
        next.push_back({-1, -1});
        cnt.push_back(0);
        return next.size() - 1;
    }

    static void check_range(unsigned long long x){
        /// BITS % 64 keeps the shift count legal in the BITS == 64 instantiation, where the left side short-circuits
        assert(BITS == 64 || x >> (BITS % 64) == 0);
    }

    long long weight(int v) const{
        return v == -1 ? 0 : cnt[v];
    }
};

int main(){
    /// Word and prefix counts
    Trie<> trie;
    for (string s : {"apple", "app", "apply", "banana", "app"}) trie.insert(s);
    assert(trie.count_prefix("app") == 4);  /// app twice, apple, apply
    assert(trie.count_word("app") == 2);
    assert(trie.count_word("appl") == 0);   /// a prefix, never inserted whole

    /// Walk one character at a time, reading each node's counts
    int v = trie.child(0, 'a');
    v = trie.child(v, 'p');
    v = trie.child(v, 'p');
    assert(v == trie.find("app"));
    assert(trie.pass[v] == 4);
    assert(trie.ends[v] == 2);
    assert(trie.child(v, 'x') == -1);

    /// Xor queries over {3, 10, 5, 25}
    BinaryTrie<5> bt;
    for (int x : {3, 10, 5, 25}) bt.insert(x);
    assert(bt.min_xor(6) == 3);             /// 6 ^ {3, 10, 5, 25} = {5, 12, 3, 31}
    assert(bt.max_xor(6) == 31);
    assert(bt.kth_xor(6, 1) == 5);

    bt.erase(5);
    assert(bt.count(5) == 0);
    assert(bt.min_xor(6) == 5);

    return 0;
}
