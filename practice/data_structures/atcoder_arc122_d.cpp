/***
 *
 * Trie and Binary Trie
 * Prefix tree that counts inserted words by prefix, and a binary trie answering xor queries over a multiset of integers
 *
 * Complexity:
 *   - Trie: O(|s|) per operation, O(total length * SIGMA) memory
 *   - BinaryTrie: O(BITS) per operation, O(distinct inserted values * BITS) memory
 *
 * Trie<SIGMA = 26, BASE = 'a'> trie: no words yet, characters must lie in [BASE, BASE + SIGMA)
 * Trie<> trie(max_nodes): the same with memory for max_nodes >= 1 nodes reserved up front, more nodes still fit
 * trie.insert(s): adds one copy of s, duplicates count separately
 * trie.count_prefix(p): inserted words that start with p, every word for the empty prefix
 * trie.count_word(s): copies of exactly s
 * trie.find(s): node of the prefix s, -1 if no inserted word starts with s, the root 0 for the empty prefix
 * trie.child(v, ch): node one character ch below node v, -1 if there is none or ch is outside the alphabet, v != -1
 * trie.pass[v]: inserted words that start with node v's prefix
 * trie.ends[v]: copies of exactly node v's prefix
 *
 * A trie holds at most 1 + total inserted length nodes, the natural max_nodes
 * Each node takes 4 * SIGMA + 8 bytes
 * Without max_nodes the vector doubling briefly holds up to twice the memory of the nodes in use
 *
 * BinaryTrie<BITS = 30> bt: an empty multiset of unsigned long long values in [0, 2^BITS), 1 <= BITS <= 64
 * bt.insert(x, c = 1): adds c >= 1 copies of x
 * bt.erase(x, c = 1): removes c >= 1 copies of x, count(x) >= c
 * bt.count(x): copies of x
 * bt.size(): copies of all values
 * bt.min_xor(x): min of x ^ y over the y in the multiset, size() > 0
 * bt.max_xor(x): max of x ^ y over the y in the multiset, size() > 0
 * bt.kth_xor(x, k): k-th smallest x ^ y counting copies, 0-based, 0 <= k < size()
 *
 * Erased nodes stay allocated, so memory grows with insertions of new values, not with live size
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

    /// Reserving up front skips the vector doubling, whose copy otherwise holds up to twice the nodes in memory at once
    explicit Trie(int max_nodes){
        assert(max_nodes >= 1);
        next.reserve(max_nodes), pass.reserve(max_nodes), ends.reserve(max_nodes);
        new_node();
    }

    int child(int v, char ch) const{
        int c = ch - BASE;
        return c < 0 || c >= SIGMA ? -1 : next[v][c];
    }

    /// Deleted so that child(v, s[i] - 'a') fails to compile instead of reading the int as a character
    int child(int v, int c) const = delete;

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

    int find(const string& s) const{
        int v = 0;
        for (int i = 0; i < (int)s.size() && v != -1; i++) v = child(v, s[i]);
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

template <int BITS = 30>
struct BinaryTrie{
    static_assert(1 <= BITS && BITS <= 64, "BITS must lie in [1, 64]");

    vector<array<int, 2>> next;
    vector<long long> cnt;

    BinaryTrie(){
        new_node();
    }

    int new_node(){
        next.push_back({-1, -1});
        cnt.push_back(0);
        return next.size() - 1;
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

    void erase(unsigned long long x, long long c = 1){
        assert(c >= 1 && count(x) >= c);
        int v = 0;
        cnt[0] -= c;
        for (int b = BITS - 1; b >= 0; b--){
            v = next[v][x >> b & 1];
            cnt[v] -= c;
        }
    }

    long long count(unsigned long long x) const{
        check_range(x);
        int v = 0;
        for (int b = BITS - 1; b >= 0 && v != -1; b--) v = next[v][x >> b & 1];
        return weight(v);
    }

    long long size() const{
        return cnt[0];
    }

    unsigned long long min_xor(unsigned long long x) const{
        return kth_xor(x, 0);
    }

    unsigned long long max_xor(unsigned long long x) const{
        return kth_xor(x, size() - 1);
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

    long long weight(int v) const{
        return v == -1 ? 0 : cnt[v];
    }

    static void check_range(unsigned long long x){
        /// BITS % 64 keeps the shift count legal in the BITS == 64 instantiation, where the left side short-circuits
        assert(BITS == 64 || x >> (BITS % 64) == 0);
    }
};

int solve(vector<int>& a, int bit){
    if (bit < 0 || a.empty()) return 0;

    vector<int> lo, hi;
    for (int x : a) (x >> bit & 1 ? hi : lo).push_back(x);
    if (lo.size() % 2 == 0) return max(solve(lo, bit - 1), solve(hi, bit - 1));

    BinaryTrie<30> bt;
    for (int x : lo) bt.insert(x);
    unsigned long long best = ULLONG_MAX;
    for (int x : hi) best = min(best, bt.min_xor(x));
    return best;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<int> a(2 * n);
    for (int& x : a) cin >> x;
    cout << solve(a, 29) << '\n';
    return 0;
}
