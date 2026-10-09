/***
 *
 * Suffix Automaton
 * Smallest DFA accepting every suffix of a string, built online one character at a time
 *
 * Complexity: O(n * SIGMA) time and memory to build, at most max(n + 1, 2n - 1) states
 * contains / occurrences: O(|pattern|), plus O(n) on the first occurrences call after an add
 * longest_common_substring(a, b): O(|a| * SIGMA + |b|)
 *
 * SuffixAutomaton<SIGMA, BASE> sa(s); characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * sa.add(c): appends c to the string
 * sa.distinct_substrings(): number of distinct non-empty substrings so far
 * sa.contains(p): true if p is a substring, the empty string included
 * sa.occurrences(p): how many positions p occurs at (overlapping), occurrences("") = n + 1
 * State v: len[v] = longest string in v, link[v] = suffix link, next[v][c] = transition or -1
 * State 0 is the root (empty string), link[0] = -1
 *
 * longest_common_substring<SIGMA, BASE>(a, b): one longest common substring, empty if none
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int SIGMA = 26, char BASE = 'a'>
struct SuffixAutomaton{
    vector<array<int, SIGMA>> next;
    vector<int> len, link, is_prefix, endpos;
    int last = 0;
    long long distinct = 0;

    SuffixAutomaton(){
        new_state(0, -1, 0);
    }

    explicit SuffixAutomaton(const string& s) : SuffixAutomaton(){
        for (char c : s) add(c);
    }

    int new_state(int length, int suffix, int prefix){
        next.push_back({});
        next.back().fill(-1);
        len.push_back(length), link.push_back(suffix), is_prefix.push_back(prefix);
        return len.size() - 1;
    }

    void add(char ch){
        int c = ch - BASE;
        assert(0 <= c && c < SIGMA);
        int cur = new_state(len[last] + 1, 0, 1), p = last;
        last = cur;

        for (; p != -1 && next[p][c] == -1; p = link[p]) next[p][c] = cur;
        if (p != -1){
            int q = next[p][c];
            if (len[p] + 1 == len[q]) link[cur] = q;
            else{
                int clone = new_state(len[p] + 1, link[q], 0);
                next[clone] = next[q];
                for (; p != -1 && next[p][c] == q; p = link[p]) next[p][c] = clone;
                link[q] = link[cur] = clone;
            }
        }

        /// The new substrings are exactly the suffixes of the string longer than the suffix link's
        distinct += len[cur] - len[link[cur]];
    }

    long long distinct_substrings() const{
        return distinct;
    }

    /// State reached by reading p from the root, -1 if p is not a substring
    int walk(const string& p) const{
        int v = 0;
        for (char ch : p){
            int c = ch - BASE;
            if (c < 0 || c >= SIGMA || next[v][c] == -1) return -1;
            v = next[v][c];
        }

        return v;
    }

    bool contains(const string& p) const{
        return walk(p) != -1;
    }

    long long occurrences(const string& p){
        if (p.empty()) return len[last] + 1;
        int v = walk(p);
        if (v == -1) return 0;

        if (endpos.size() != len.size()) build_endpos();
        return endpos[v];
    }

    /// |endpos(v)| = prefixes ending in v's subtree of the suffix link tree, links point to shorter states
    void build_endpos(){
        int states = len.size(), n = len[last];
        vector<int> bucket(n + 2, 0), order(states);
        for (int v = 0; v < states; v++) bucket[len[v] + 1]++;
        for (int i = 1; i <= n + 1; i++) bucket[i] += bucket[i - 1];
        for (int v = 0; v < states; v++) order[bucket[len[v]]++] = v;

        endpos.assign(is_prefix.begin(), is_prefix.end());
        for (int i = states - 1; i > 0; i--) endpos[link[order[i]]] += endpos[order[i]];
    }
};

template <int SIGMA = 26, char BASE = 'a'>
string longest_common_substring(const string& a, const string& b){
    SuffixAutomaton<SIGMA, BASE> sa(a);
    int v = 0, matched = 0, best = 0, best_end = 0;

    for (int i = 0; i < (int)b.size(); i++){
        int c = b[i] - BASE;
        if (c < 0 || c >= SIGMA){
            v = 0, matched = 0;
            continue;
        }

        while (v != 0 && sa.next[v][c] == -1) v = sa.link[v], matched = sa.len[v];
        if (sa.next[v][c] != -1) v = sa.next[v][c], matched++;
        if (matched > best) best = matched, best_end = i + 1;
    }

    return b.substr(best_end - best, best);
}

int main(){
    SuffixAutomaton<> sa("abcbc");
    assert(sa.distinct_substrings() == 12);
    assert(sa.contains("") && sa.contains("cbc") && sa.contains("abcbc") && sa.contains("bcb"));
    assert(!sa.contains("ac") && !sa.contains("abcbcb") && !sa.contains("A") && !sa.contains("c{"));
    assert(sa.occurrences("bc") == 2 && sa.occurrences("c") == 2 && sa.occurrences("a") == 1);
    assert(sa.occurrences("cbc") == 1 && sa.occurrences("ca") == 0 && sa.occurrences("") == 6);

    sa.add('b');
    sa.add('c');
    assert(sa.distinct_substrings() == 18);
    assert(sa.occurrences("bc") == 3 && sa.occurrences("bcbc") == 2 && sa.occurrences("cb") == 2);

    SuffixAutomaton<> same("aaaa");
    assert(same.distinct_substrings() == 4 && (int)same.len.size() == 5);
    assert(same.occurrences("a") == 4 && same.occurrences("aa") == 3 && same.occurrences("aaaa") == 1);

    SuffixAutomaton<> empty;
    assert(empty.distinct_substrings() == 0 && empty.contains("") && !empty.contains("a"));
    assert(empty.occurrences("") == 1 && empty.occurrences("a") == 0);

    SuffixAutomaton<> skewed("abbb");
    assert((int)skewed.len.size() == 7 && skewed.distinct_substrings() == 7);

    assert(longest_common_substring("xabcdey", "zzabcdq") == "abcd");
    assert(longest_common_substring("banana", "ananas") == "anana");
    assert(longest_common_substring("abc", "def") == "");
    assert(longest_common_substring("", "abc") == "");
    assert(longest_common_substring("abc", "") == "");
    assert(longest_common_substring("abab", "b?aba") == "aba");

    SuffixAutomaton<2, '0'> binary("0110");
    assert(binary.distinct_substrings() == 8 && binary.occurrences("1") == 2 && !binary.contains("00"));

    return 0;
}
