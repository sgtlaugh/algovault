/***
 *
 * Subsequence Automaton
 * Next occurrence of a symbol at or after a position, for greedy subsequence matching
 *
 * Complexity:
 *   - SubsequenceAutomaton: O(n alphabet) build time and memory, O(1) per next, O(|t|) per is_subsequence
 *   - SparseSubsequenceAutomaton: O(n log n) build, O(log n) per next, O(|t| log n) per is_subsequence
 *   - distinct_subsequences: O(n log n)
 *
 * SubsequenceAutomaton sa(s, alphabet = 26, first = 'a'): symbols of s must lie in [first, first + alphabet)
 * sa.is_subsequence(t): true if t is a subsequence of s, the empty t included
 * sa.next(i, c): smallest j >= i with s[j] = c, n if none, for 0 <= i <= n + 1
 *
 * SparseSubsequenceAutomaton sa(s): the same calls for a vector s of any values, such as ints up to 1e9
 *
 * distinct_subsequences(s, mod = 998244353): number of distinct non-empty subsequences of s, modulo mod
 *
 * Notes:
 *   - next(n + 1, c) = n, so a greedy walk i = next(i, c) + 1 never leaves the table and matched t iff i <= n
 *   - A symbol c outside the alphabet never occurs, so next returns n
 *   - The dense table holds 4 (n + 2) alphabet bytes, 104 MB at n = 1e6 and alphabet 26; above that use the sparse one
 *   - s is a string, a C string or a vector; distinct_subsequences needs mod <= 2^62
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SubsequenceAutomaton{
    int n, alphabet, first;
    vector<int> next_pos;

    template<typename Container>
    SubsequenceAutomaton(const Container& s, int alphabet = 26, int first = 'a')
        : n(s.size()), alphabet(alphabet), first(first), next_pos((n + 2) * alphabet, n){
        for (int i = n - 1; i >= 0; i--){
            long long c = s[i] - (long long)first;
            assert(0 <= c && c < alphabet);
            copy_n(next_pos.begin() + (i + 1) * alphabet, alphabet, next_pos.begin() + i * alphabet);
            next_pos[i * alphabet + c] = i;
        }
    }

    SubsequenceAutomaton(const char* s, int alphabet = 26, int first = 'a')
        : SubsequenceAutomaton(string(s), alphabet, first) {}

    template<typename Container>
    bool is_subsequence(const Container& t) const{
        int i = 0;

        for (long long c : t){
            if (i > n) return false;
            i = next(i, c) + 1;
        }

        return i <= n;
    }

    bool is_subsequence(const char* t) const{
        return is_subsequence(string(t));
    }

    /// long long so a large value is never narrowed into the alphabet
    int next(int i, long long c) const{
        c -= first;
        if (c < 0 || c >= alphabet) return n;

        return next_pos[i * alphabet + c];
    }
};

template<typename T>
struct SparseSubsequenceAutomaton{
    int n;
    vector<pair<T, int>> occurrences;

    SparseSubsequenceAutomaton(const vector<T>& s) : n(s.size()){
        for (int i = 0; i < n; i++) occurrences.emplace_back(s[i], i);
        sort(occurrences.begin(), occurrences.end());
    }

    bool is_subsequence(const vector<T>& t) const{
        int i = 0;

        for (const T& c : t){
            if (i > n) return false;
            i = next(i, c) + 1;
        }

        return i <= n;
    }

    int next(int i, const T& c) const{
        auto it = lower_bound(occurrences.begin(), occurrences.end(), make_pair(c, i));

        return it != occurrences.end() && it->first == c ? it->second : n;
    }
};

template<typename Container>
long long distinct_subsequences(const Container& s, long long mod = 998244353){
    int n = s.size();
    vector<int> order(n), prev(n, -1);

    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){ return make_pair(s[a], a) < make_pair(s[b], b); });
    for (int k = 1; k < n; k++){
        if (s[order[k]] == s[order[k - 1]]) prev[order[k]] = order[k - 1];
    }

    /// count[i] = distinct subsequences of s[0..i), the empty one included; those ending at s[i] were new after prev
    vector<long long> count(n + 1);
    count[0] = 1;
    for (int i = 0; i < n; i++){
        long long repeated = prev[i] >= 0 ? count[prev[i]] : 0;
        count[i + 1] = ((2 * count[i] - repeated) % mod + mod) % mod;
    }

    return (count[n] - 1 + mod) % mod;
}

long long distinct_subsequences(const char* s, long long mod = 998244353){
    return distinct_subsequences(string(s), mod);
}

int main(){
    /// "banana": next jumps to the next copy of a symbol, n = 6 when none is left
    SubsequenceAutomaton sa("banana");
    assert(sa.next(0, 'a') == 1);
    assert(sa.next(2, 'a') == 3);
    assert(sa.next(4, 'b') == 6);
    assert(sa.next(0, 'z') == 6);
    assert(sa.is_subsequence("bnn"));
    assert(sa.is_subsequence(""));
    assert(!sa.is_subsequence("nb"));

    SubsequenceAutomaton bits(vector<int>{1, 0, 1}, 2, 0);
    assert(bits.is_subsequence(vector<int>{1, 1}));
    assert(!bits.is_subsequence(vector<int>{0, 0}));

    /// Large values, as Library Checker number_of_subsequences gives them
    SparseSubsequenceAutomaton big(vector<int>{1000000000, 7, 1000000000});
    assert(big.next(1, 1000000000) == 2);
    assert(big.next(0, 8) == 3);
    assert(big.is_subsequence(vector<int>{7, 1000000000}));
    assert(!big.is_subsequence(vector<int>{7, 7}));

    assert(distinct_subsequences("aba") == 6);  /// a, b, ab, ba, aa, aba
    assert(distinct_subsequences(vector<int>{0, 0, 0, 1, 1}) == 11);
    assert(distinct_subsequences("") == 0);

    return 0;
}
