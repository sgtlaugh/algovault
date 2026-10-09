/***
 *
 * Dynamic String Hash
 * Polynomial hash of any substring while ranges of the string are reassigned to a single character
 *
 * Complexity: O(n) to build, O(log n) per assign and per hash query
 *
 * DynamicStringHash h(s): 0-based positions, ranges [l, r] inclusive
 * h.assign(l, r, c): sets every character in [l, r] to c
 * h.hash(l, r): hash of s[l..r] modulo 2^61 - 1, equal substrings give equal hashes
 *
 * The base is random per run (as in hashing.cpp) so anti-hash tests cannot target it, pass a seed to fix it
 * Instances built with the default seed share a base, so their hashes are comparable
 * Two different substrings of equal length collide with probability about length / 2^61
 * Memory: about 64n bytes (segment tree and lazy tags over 4n nodes, plus power tables)
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicStringHash{
    static constexpr unsigned long long MOD = (1ULL << 61) - 1;

    int n;
    unsigned long long base;
    vector<unsigned long long> power, run;  /// run[k] = 1 + base + ... + base^(k - 1)
    vector<unsigned long long> tree;
    vector<int> lazy;                       /// pending character + 1, 0 for none

    static unsigned long long mul(unsigned long long a, unsigned long long b){
        unsigned __int128 x = (unsigned __int128)a * b;
        unsigned long long res = (unsigned long long)(x & MOD) + (unsigned long long)(x >> 61);
        return res >= MOD ? res - MOD : res;
    }

    static unsigned long long add(unsigned long long a, unsigned long long b){
        a += b;
        return a >= MOD ? a - MOD : a;
    }

    static inline const unsigned long long default_seed = chrono::steady_clock::now().time_since_epoch().count();

    static unsigned long long value(int c){
        return c + 1;
    }

    DynamicStringHash(const string& s, unsigned long long seed = default_seed)
        : n(s.size()), power(s.size() + 1), run(s.size() + 1), tree(4 * max<size_t>(1, s.size())), lazy(4 * max<size_t>(1, s.size()), 0){
        base = mt19937_64(seed)() % (MOD - 1000) + 500;
        power[0] = 1, run[0] = 0;
        for (int i = 1; i <= n; i++) power[i] = mul(power[i - 1], base), run[i] = add(mul(run[i - 1], base), 1);
        if (n) build(1, 0, n - 1, s);
    }

    void build(int node, int a, int b, const string& s){
        if (a == b){
            tree[node] = value((unsigned char)s[a]);
            return;
        }

        int m = (a + b) / 2;
        build(2 * node, a, m, s);
        build(2 * node + 1, m + 1, b, s);
        tree[node] = add(mul(tree[2 * node], power[b - m]), tree[2 * node + 1]);
    }

    void apply(int node, int len, int c){
        tree[node] = mul(value(c), run[len]);
        lazy[node] = c + 1;
    }

    void push(int node, int a, int b){
        if (!lazy[node]) return;
        int m = (a + b) / 2;
        apply(2 * node, m - a + 1, lazy[node] - 1);
        apply(2 * node + 1, b - m, lazy[node] - 1);
        lazy[node] = 0;
    }

    void assign(int node, int a, int b, int l, int r, int c){
        if (r < a || b < l) return;
        if (l <= a && b <= r) return apply(node, b - a + 1, c);
        push(node, a, b);
        int m = (a + b) / 2;
        assign(2 * node, a, m, l, r, c);
        assign(2 * node + 1, m + 1, b, l, r, c);
        tree[node] = add(mul(tree[2 * node], power[b - m]), tree[2 * node + 1]);
    }

    /// Hash of [max(a, l), min(b, r)] inside this node
    unsigned long long query(int node, int a, int b, int l, int r){
        if (r < a || b < l) return 0;
        if (l <= a && b <= r) return tree[node];
        push(node, a, b);
        int m = (a + b) / 2;
        unsigned long long left = query(2 * node, a, m, l, r), right = query(2 * node + 1, m + 1, b, l, r);
        int right_len = max(0, min(b, r) - m);  /// when the left part is non-empty the right part starts at m + 1
        return add(mul(left, power[right_len]), right);
    }

    void assign(int l, int r, char c){
        assert(0 <= l && l <= r && r < n);
        assign(1, 0, n - 1, l, r, (unsigned char)c);
    }

    unsigned long long hash(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r);
    }
};

int main(){
    DynamicStringHash h("abcabcab");
    assert(h.hash(0, 2) == h.hash(3, 5));
    assert(h.hash(0, 1) == h.hash(6, 7));
    assert(h.hash(0, 2) != h.hash(1, 3));

    h.assign(3, 5, 'z');
    assert(h.hash(0, 2) != h.hash(3, 5));
    assert(h.hash(3, 4) == h.hash(4, 5));

    DynamicStringHash p("abcabcab", 11), q("abczzzab", 11);
    p.assign(3, 5, 'z');
    assert(p.hash(0, 7) == q.hash(0, 7));

    DynamicStringHash a("xyxy", 7), b("yxyx", 7);
    assert(a.hash(0, 1) == b.hash(1, 2) && a.hash(1, 3) == b.hash(0, 2));
    a.assign(0, 3, 'q'), b.assign(1, 2, 'q');
    assert(a.hash(1, 2) == b.hash(1, 2));

    DynamicStringHash pizza("pizza"), another_pizza("pizza");
    assert(pizza.hash(0, 4) == another_pizza.hash(0, 4));

    DynamicStringHash single("k");
    single.assign(0, 0, 'm');
    assert(single.hash(0, 0) == (unsigned long long)'m' + 1);

    return 0;
}
