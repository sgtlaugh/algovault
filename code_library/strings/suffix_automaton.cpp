/***
 *
 * Suffix Automaton
 * Smallest DFA accepting every substring of one string or of several, built online one character at a time
 *
 * Complexity, for total length n over all strings:
 *   - build: O(n SIGMA) time and memory, at most max(n + 1, 2n - 1) states
 *   - contains, occurrences: O(|pattern|), plus O(n) on the first occurrences call after an add
 *   - string_counts: O(n sqrt(n)) worst case, O(n) for one string
 *   - longest_common_substring(a, b): O(|a| SIGMA + |b|)
 *
 * SuffixAutomaton<SIGMA = 26, BASE = 'a'> sa(s): the automaton of s, characters must lie in [BASE, BASE + SIGMA)
 * sa.add(c): appends c to the current string
 * sa.add_string(s): starts a new string holding s, so substrings never span two strings, nothing for an empty s
 * sa.contains(p): true if p is a substring of some string, the empty string included
 * sa.distinct_substrings(): number of distinct non-empty substrings over all strings so far
 * sa.occurrences(p): how many positions p occurs at over all strings (overlapping), occurrences("") = n + strings
 * sa.string_counts(strs): res[v] = how many of strs contain the strings of state v, strs exactly the added strings
 * sa.strings: number of non-empty strings, at least 1
 * State v: len[v] = longest string in v, link[v] = suffix link, next[v][c] = transition or -1
 * State 0 is the root (empty string), link[0] = -1
 *
 * longest_common_substring<SIGMA, BASE>(a, b): one longest common substring, empty if none
 *
 * Notes:
 *   - Strings shared by k of m strings: states with string_counts >= k, each adding len[v] - len[link[v]] strings
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int SIGMA = 26, char BASE = 'a'>
struct SuffixAutomaton{
    vector<array<int, SIGMA>> next;
    vector<int> len, link, prefix_ends, endpos;
    int last = 0, strings = 1;
    long long distinct = 0, chars = 0, endpos_chars = -1;

    SuffixAutomaton(){
        new_state(0, -1, 0);
    }

    explicit SuffixAutomaton(const string& s) : SuffixAutomaton(){
        size_t states = 2 * s.size() + 1;
        next.reserve(states), len.reserve(states), link.reserve(states), prefix_ends.reserve(states);

        for (char c : s) add(c);
    }

    void add(char ch){
        int c = ch - BASE;

        assert(0 <= c && c < SIGMA);
        chars++;

        /// Only after add_string: the new prefix is already a substring, so it gets no new state, at most a clone
        if (next[last][c] != -1){
            int q = next[last][c];
            last = len[last] + 1 == len[q] ? q : split(last, c);
            prefix_ends[last]++;
            return;
        }

        int cur = new_state(len[last] + 1, 0, 1), p = last;
        last = cur;

        for (; p != -1 && next[p][c] == -1; p = link[p]) next[p][c] = cur;
        if (p != -1) link[cur] = len[p] + 1 == len[next[p][c]] ? next[p][c] : split(p, c);

        /// The new substrings are exactly the suffixes of the string longer than the suffix link's
        distinct += len[cur] - len[link[cur]];
    }

    void add_string(const string& s){
        if (s.empty()) return;
        if (chars > 0) strings++;
        last = 0;
        for (char c : s) add(c);
    }

    bool contains(const string& p) const{
        return walk(p) != -1;
    }

    long long distinct_substrings() const{
        return distinct;
    }

    long long occurrences(const string& p){
        if (p.empty()) return chars + strings;
        int v = walk(p);
        if (v == -1) return 0;

        if (endpos_chars != chars) build_endpos();
        return endpos[v];
    }

    /// Each string's prefix states climb their suffix links once per string, O(L sqrt(L)) for total length L
    vector<int> string_counts(const vector<string>& strs) const{
        vector<int> res(len.size(), 0), seen(len.size(), -1);
        for (int id = 0; id < (int)strs.size(); id++){
            int v = 0;
            for (char ch : strs[id]){
                int c = ch - BASE;
                assert(0 <= c && c < SIGMA);
                v = next[v][c];
                assert(v != -1);
                for (int u = v; u != -1 && seen[u] != id; u = link[u]) seen[u] = id, res[u]++;
            }
        }

        return res;
    }

    int new_state(int length, int suffix, int prefix){
        next.push_back({});
        next.back().fill(-1);
        len.push_back(length), link.push_back(suffix), prefix_ends.push_back(prefix);
        return len.size() - 1;
    }

    /// Gives q = next[p][c] a copy of length len[p] + 1 that takes over q's shorter strings, returns the copy
    int split(int p, int c){
        int q = next[p][c], clone = new_state(len[p] + 1, link[q], 0);
        next[clone] = next[q];
        for (; p != -1 && next[p][c] == q; p = link[p]) next[p][c] = clone;
        link[q] = clone;
        return clone;
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

    /// |endpos(v)| = prefix ends in v's subtree of the suffix link tree, links point to shorter states
    void build_endpos(){
        int states = len.size(), n = *max_element(len.begin(), len.end());
        vector<int> bucket(n + 2, 0), order(states);

        for (int v = 0; v < states; v++) bucket[len[v] + 1]++;
        for (int i = 1; i <= n + 1; i++) bucket[i] += bucket[i - 1];
        for (int v = 0; v < states; v++) order[bucket[len[v]]++] = v;

        endpos.assign(prefix_ends.begin(), prefix_ends.end());
        for (int i = states - 1; i > 0; i--) endpos[link[order[i]]] += endpos[order[i]];
        endpos_chars = chars;
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
    assert(sa.distinct_substrings() == 12);  /// 15 substrings minus the repeats b, c, bc
    assert(sa.contains("cbc"));
    assert(!sa.contains("ac"));
    assert(sa.occurrences("bc") == 2);

    sa.add('b');                             /// online: now "abcbcb"
    assert(sa.occurrences("cb") == 2);

    assert(longest_common_substring("banana", "ananas") == "anana");

    /// "ab" and "bab": substrings a, b, ab, ba, bab, of which a, b and ab lie in both
    SuffixAutomaton<> both("ab");
    both.add_string("bab");
    assert(both.strings == 2);
    assert(both.distinct_substrings() == 5);
    assert(!both.contains("abb"));
    assert(both.occurrences("b") == 3);
    vector<int> in_both = both.string_counts({"ab", "bab"});
    int shared = 0;
    for (int v = 1; v < (int)in_both.size(); v++) if (in_both[v] == 2) shared += both.len[v] - both.len[both.link[v]];
    assert(shared == 3);
    return 0;
}
