/***
 *
 * 64-bit hashing for vectors or strings
 * Get the forward and reverse hash of any segment
 * Base is chosen randomly to prevent anti-hash cases from being constructed
 *
 * Complexity - O(n) to build, O(1) for each hash query, O(n) memory per instance
 *
 * PolyHash h(s): s is a string, vector or const char*, 0-based positions, ranges [l, r] inclusive
 * h.get_hash(l, r): hash of s[l..r] modulo 2^61 - 1, equal segments give equal hashes
 * h.rev_hash(l, r): hash of s[l..r] read backwards, s[l..r] is a palindrome iff get_hash(l, r) == rev_hash(l, r)
 * All instances share one base per run, so segments of different strings can be compared
 * Two different segments of equal length collide with probability about length / 2^61
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PolyHash{
    static constexpr uint64_t mod = (1ULL << 61) - 1;
    static inline const uint64_t base = mt19937_64(chrono::system_clock::now().time_since_epoch().count())() % (mod / 3) + (mod / 3);

    /// Remove suff vector and usage if reverse hash is not required for more speed
    vector<int64_t> pref, suff, base_pow;

    PolyHash() {}

    template <typename Container>
    PolyHash(const Container& ar){
        int n = ar.size();
        pref.resize(n + 3, 0), suff.resize(n + 3, 0), base_pow.resize(n + 1, 1);

        for (int i = 1; i <= n; i++){
            base_pow[i] = modmul(base_pow[i - 1], base);
            pref[i] = modmul(pref[i - 1], base) + value_of(ar[i - 1]);
            if (pref[i] >= (int64_t)mod) pref[i] -= mod;
        }

        for (int i = n; i >= 1; i--){
            suff[i] = modmul(suff[i + 1], base) + value_of(ar[i - 1]);
            if (suff[i] >= (int64_t)mod) suff[i] -= mod;
        }
    }

    PolyHash(const char* str) : PolyHash(string(str)) {}
    PolyHash(char* str) : PolyHash((const char*)str) {}

    /// a, b < mod
    static int64_t modmul(uint64_t a, uint64_t b){
        __uint128_t c = (__uint128_t)a * b;
        uint64_t r = (c & mod) + (c >> 61);
        return r >= mod ? r - mod : r;
    }

    /// The 997 offset keeps zeros from hashing like an empty prefix
    template <typename T>
    static uint64_t value_of(const T& x){
        return (uint64_t)((((__int128)x % (__int128)mod) + mod + 997) % mod);
    }

    uint64_t get_hash(int l, int r) const{
        int64_t h = pref[r + 1] - modmul(base_pow[r - l + 1], pref[l]);
        return h < 0 ? h + mod : h;
    }

    uint64_t rev_hash(int l, int r) const{
        int64_t h = suff[l + 1] - modmul(base_pow[r - l + 1], suff[r + 2]);
        return h < 0 ? h + mod : h;
    }
};

/***
 *
 * Palindromic Tree (Eertree)
 * One node per distinct palindromic substring, built online one character at a time
 *
 * Complexity: O(n * SIGMA) memory, amortized O(1) per add, O(n log n) for min_palindromic_factorization
 *
 * PalindromicTree<SIGMA, BASE> tree; characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * tree.add(c): appends c, returns true if a new distinct palindrome appeared
 * tree.distinct(): number of distinct non-empty palindromic substrings so far
 * tree.longest: length of the longest palindromic substring so far
 * tree.occurrences(): occurrences[v] = how many times node v occurs as a substring, nodes 2.. are real palindromes
 * Node v: len[v], link[v] = longest proper palindromic suffix, next[v][c] = c + palindrome(v) + c
 * Nodes 0 and 1 are the roots: node 0 is the empty palindrome, node 1 the imaginary one of length -1
 *
 * min_palindromic_factorization(s): res[i] = fewest palindromes that the prefix s[0..i] splits into
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int SIGMA = 26, char BASE = 'a'>
struct PalindromicTree{
    string s;
    vector<array<int, SIGMA>> next;
    vector<int> len, link, diff, series_link, ends_here;
    int last = 0, longest = 0;

    PalindromicTree(){
        new_node(0, 1);
        new_node(-1, 0);
    }

    int new_node(int length, int suffix){
        next.push_back({});
        next.back().fill(-1);
        len.push_back(length), link.push_back(suffix), ends_here.push_back(0);
        int v = len.size() - 1;
        diff.push_back(length - (v < 2 ? length : len[suffix]));
        series_link.push_back(v < 2 || diff[v] != diff[suffix] ? suffix : series_link[suffix]);
        return v;
    }

    int extend_from(int v) const{
        int i = s.size() - 1;
        while (i - len[v] - 1 < 0 || s[i - len[v] - 1] != s[i]) v = link[v];
        return v;
    }

    bool add(char ch){
        int c = ch - BASE;
        assert(0 <= c && c < SIGMA);
        s.push_back(ch);

        int v = extend_from(last);
        bool fresh = next[v][c] == -1;
        if (fresh){
            int suffix = len[v] == -1 ? 0 : next[extend_from(link[v])][c];
            int u = new_node(len[v] + 2, suffix);
            next[v][c] = u;
        }
        last = next[v][c];
        ends_here[last]++;
        longest = max(longest, len[last]);
        return fresh;
    }

    int distinct() const{
        return len.size() - 2;
    }

    /// Suffix links always point to older nodes, so one backwards pass pushes counts up
    vector<long long> occurrences() const{
        vector<long long> res(ends_here.begin(), ends_here.end());
        for (int v = len.size() - 1; v >= 2; v--) res[link[v]] += res[v];
        return res;
    }
};

template <int SIGMA = 26, char BASE = 'a'>
vector<int> min_palindromic_factorization(const string& s){
    int n = s.size();
    PalindromicTree<SIGMA, BASE> tree;
    vector<int> best(n + 1, INT_MAX), series(2, 0), res(n);
    best[0] = 0;

    for (int i = 1; i <= n; i++){
        tree.add(s[i - 1]);
        series.resize(tree.len.size(), 0);
        /// Palindromic suffixes split into O(log n) runs with equal length differences, each run reuses its link's answer
        for (int v = tree.last; tree.len[v] > 0; v = tree.series_link[v]){
            series[v] = best[i - (tree.len[tree.series_link[v]] + tree.diff[v])];
            if (tree.diff[v] == tree.diff[tree.link[v]]) series[v] = min(series[v], series[tree.link[v]]);
            best[i] = min(best[i], series[v] + 1);
        }
        res[i - 1] = best[i];
    }

    return res;
}

#include <malloc.h>

int main(){
    /// Each case allocates and frees ~60 MB of tree nodes; keeping freed memory in the heap reuses those pages instead of
    /// returning them to the kernel and faulting fresh zeroed ones in for the next case (halves the run time on 40 max cases)
    mallopt(M_MMAP_MAX, 0), mallopt(M_TRIM_THRESHOLD, -1);

    int t;
    scanf("%d", &t);
    static char word[300005];
    while (t--){
        scanf("%s", word);
        int n = strlen(word);
        PolyHash hash(word);

        /// Letters packed into 0..51 keep the tree at 52 children per node instead of the 58 from 'A' to 'z'
        PalindromicTree<52, 0> tree;
        tree.next.reserve(n + 2);

        int res = 0;
        for (int i = 0; i < n; i++){
            char c = islower(word[i]) ? word[i] - 'a' : word[i] - 'A' + 26;
            if (!tree.add(c)) continue;

            /// A palindrome of the form XX is w wR w wR, as X = reverse(X) follows from XX being a palindrome
            int len = tree.len[tree.last], half = len / 2, l = i - len + 1;
            if (len % 4 == 0 && hash.get_hash(l, l + half - 1) == hash.get_hash(l + half, i)) res = max(res, len);
        }
        printf("%d\n", res);
    }
    return 0;
}
