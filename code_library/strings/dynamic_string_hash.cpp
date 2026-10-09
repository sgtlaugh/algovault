/***
 *
 * Dynamic String Hash
 * Forward and reverse polynomial hash of any substring while ranges of the string are reassigned to a single character
 *
 * Complexity: O(n) to build, O(log n) per assign and per hash query
 *
 * DynamicStringHash h(s): 0-based positions, ranges [l, r] inclusive
 * h.assign(l, r, c): sets every character in [l, r] to c, h.assign(i, i, c) is a point update
 * h.hash(l, r): hash of s[l..r] modulo 2^61 - 1, equal substrings give equal hashes
 * h.rev_hash(l, r): hash of s[l..r] read backwards, equal to hash() of the same characters in reverse order
 * h.is_palindrome(l, r): whether s[l..r] reads the same backwards, one query
 *
 * The base is random per run (as in hashing.cpp) so anti-hash tests cannot target it, pass a seed to fix it
 * Instances built with the default seed share a base, so their hashes are comparable
 * Two different substrings of equal length collide with probability about length / 2^61
 * Memory: about 96n bytes (two hash trees and lazy tags over 4n nodes, plus power tables)
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicStringHash{
    static constexpr unsigned long long MOD = (1ULL << 61) - 1;

    struct Part{
        unsigned long long fwd, rev;
        int len;
    };

    int n;
    unsigned long long base;
    vector<unsigned long long> power, run;  /// run[k] = 1 + base + ... + base^(k - 1)
    vector<unsigned long long> fwd, rev;    /// hash of each node's range read forwards and backwards
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
        : n(s.size()), power(s.size() + 1), run(s.size() + 1), fwd(4 * max<size_t>(1, s.size())),
          rev(4 * max<size_t>(1, s.size())), lazy(4 * max<size_t>(1, s.size()), 0){
        base = mt19937_64(seed)() % (MOD - 1000) + 500;
        power[0] = 1, run[0] = 0;
        for (int i = 1; i <= n; i++) power[i] = mul(power[i - 1], base), run[i] = add(mul(run[i - 1], base), 1);
        if (n) build(1, 0, n - 1, s);
    }

    /// The reversed concatenation is the reversed right part followed by the reversed left part
    Part combine(const Part& x, const Part& y) const{
        return {add(mul(x.fwd, power[y.len]), y.fwd), add(mul(y.rev, power[x.len]), x.rev), x.len + y.len};
    }

    void pull(int node, int left_len, int right_len){
        fwd[node] = add(mul(fwd[2 * node], power[right_len]), fwd[2 * node + 1]);
        rev[node] = add(mul(rev[2 * node + 1], power[left_len]), rev[2 * node]);
    }

    void build(int node, int a, int b, const string& s){
        if (a == b){
            fwd[node] = rev[node] = value((unsigned char)s[a]);
            return;
        }

        int m = (a + b) / 2;
        build(2 * node, a, m, s);
        build(2 * node + 1, m + 1, b, s);
        pull(node, m - a + 1, b - m);
    }

    void apply(int node, int len, int c){
        fwd[node] = rev[node] = mul(value(c), run[len]);
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
        pull(node, m - a + 1, b - m);
    }

    Part query(int node, int a, int b, int l, int r){
        if (r < a || b < l) return {0, 0, 0};
        if (l <= a && b <= r) return {fwd[node], rev[node], b - a + 1};
        push(node, a, b);
        int m = (a + b) / 2;
        return combine(query(2 * node, a, m, l, r), query(2 * node + 1, m + 1, b, l, r));
    }

    Part query(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r);
    }

    void assign(int l, int r, char c){
        assert(0 <= l && l <= r && r < n);
        assign(1, 0, n - 1, l, r, (unsigned char)c);
    }

    unsigned long long hash(int l, int r){
        return query(l, r).fwd;
    }

    unsigned long long rev_hash(int l, int r){
        return query(l, r).rev;
    }

    bool is_palindrome(int l, int r){
        Part p = query(l, r);
        return p.fwd == p.rev;
    }
};

int main(){
    DynamicStringHash h("abcabcab");
    assert(h.hash(0, 2) == h.hash(3, 5));  /// "abc" twice
    assert(h.hash(0, 2) != h.hash(1, 3));  /// "abc" against "bca"

    h.assign(3, 5, 'z');                   /// abczzzab
    assert(h.hash(0, 2) != h.hash(3, 5));
    assert(h.hash(3, 4) == h.hash(4, 5));  /// "zz" twice

    DynamicStringHash word("racecat");
    assert(word.is_palindrome(1, 5));                /// "aceca"
    assert(!word.is_palindrome(0, 6));
    assert(word.hash(0, 2) != word.rev_hash(4, 6));  /// "rac" against "cat" read backwards
    word.assign(6, 6, 'r');                          /// racecar
    assert(word.is_palindrome(0, 6));
    assert(word.hash(0, 2) == word.rev_hash(4, 6));  /// "rac" against "car" read backwards
    return 0;
}
