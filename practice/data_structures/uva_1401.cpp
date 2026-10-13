/***
 *
 * Trie and Binary Trie
 * Trie: prefix tree over strings, counting how many inserted words pass through or end at each node
 * BinaryTrie: multiset of BITS-bit unsigned integers answering xor queries against it
 *
 * Complexity: Trie O(|s|) per operation, O(total length * SIGMA) memory
 *             BinaryTrie O(BITS) per operation, O(distinct inserted values * BITS) memory
 *
 * Trie<SIGMA, BASE> trie; characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * trie.insert(s): adds one copy of s, duplicates count separately
 * trie.count_prefix(p): inserted words that start with p (every word for the empty prefix)
 * trie.count_word(s): copies of exactly s
 *
 * BinaryTrie<BITS> bt; values are unsigned long long in [0, 2^BITS), 1 <= BITS <= 64, defaults to 30
 * bt.insert(x, c = 1), bt.erase(x, c = 1): add or remove c >= 1 copies of x, erase requires count(x) >= c
 * bt.count(x): copies of x, bt.size(): copies of all values
 * bt.min_xor(x), bt.max_xor(x): min / max of x ^ y over the y in the multiset, requires size() > 0
 * bt.kth_xor(x, k): k-th smallest (0-indexed) of x ^ y counting copies, requires 0 <= k < size()
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

const int MOD = 20071027;

int main(){
    static char buf[300005], word[105];
    for (int tc = 1; scanf("%s", buf) == 1; tc++){
        int m;
        scanf("%d", &m);
        Trie<> trie;
        for (int i = 0; i < m; i++){
            scanf("%s", word);
            trie.insert(word);
        }

        /// dp[i] counts the splits of the suffix from i: walk the trie along it, each word end adds dp of what follows
        string s = buf;
        int n = s.size();
        vector<int> dp(n + 1, 0);
        dp[n] = 1;
        for (int i = n - 1; i >= 0; i--){
            for (int j = i, v = 0; j < n; j++){
                v = trie.next[v][s[j] - 'a'];
                if (v == -1) break;
                if (trie.ends[v]) dp[i] = (dp[i] + dp[j + 1]) % MOD;
            }
        }
        printf("Case %d: %d\n", tc, dp[0]);
    }
    return 0;
}
