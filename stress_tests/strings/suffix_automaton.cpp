#include "../common.h"

#define main library_main
#include "../../code_library/strings/suffix_automaton.cpp"
#undef main

string random_string(int n, int sigma){
    string s;
    for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, sigma - 1));
    return s;
}

long long count_occurrences(const string& s, const string& p){
    long long res = 0;
    for (int i = 0; i + (int)p.size() <= (int)s.size(); i++) res += s.compare(i, p.size(), p) == 0;
    return res;
}

int lcs_length(const string& a, const string& b){
    int best = 0;
    vector<vector<int>> dp(a.size() + 1, vector<int>(b.size() + 1, 0));
    for (int i = 1; i <= (int)a.size(); i++){
        for (int j = 1; j <= (int)b.size(); j++){
            if (a[i - 1] == b[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1, best = max(best, dp[i][j]);
        }
    }

    return best;
}

void check_lcs(const string& a, const string& b){
    string common = longest_common_substring<3, 'a'>(a, b);
    assert((int)common.size() == lcs_length(a, b));
    assert(a.find(common) != string::npos && b.find(common) != string::npos);
}

/// Distinct counts after every add, state bound, contains and occurrences against the set of all substrings
void check(const string& s){
    SuffixAutomaton<3, 'a'> sa;
    set<string> seen;
    for (int i = 0; i < (int)s.size(); i++){
        sa.add(s[i]);
        for (int l = 0; l <= i; l++) seen.insert(s.substr(l, i - l + 1));
        assert(sa.distinct_substrings() == (long long)seen.size());
        assert((int)sa.len.size() <= max(i + 2, 2 * (i + 1) - 1));
    }

    assert(sa.contains("") && sa.occurrences("") == (long long)s.size() + 1);
    for (const string& p : seen) assert(sa.contains(p) && sa.occurrences(p) == count_occurrences(s, p));
    for (int k = 0; k < 20; k++){
        string p = random_string(stress::rand_int(1, 6), 4);
        bool present = seen.count(p);
        assert(sa.contains(p) == present && sa.occurrences(p) == (present ? count_occurrences(s, p) : 0));
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int sigma = stress::rand_int(1, 3);
        check(random_string(stress::rand_int(0, 30), sigma));
        check_lcs(random_string(stress::rand_int(0, 30), sigma), random_string(stress::rand_int(0, 30), sigma));

        /// 'd' lies outside SIGMA = 3, so it acts as a separator that must reset the match
        check_lcs(random_string(stress::rand_int(0, 30), sigma), random_string(stress::rand_int(0, 30), 4));
    }

    for (int n = 0; n <= 10; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s;
            for (int i = 0; i < n; i++) s += mask >> i & 1 ? 'b' : 'a';
            check(s);
            check_lcs(s, "abba");
        }
    }

    /// Occurrences interleaved with adds, so the cached endpos sizes must be rebuilt
    for (long long it = 0; it < stress::scaled(300); it++){
        SuffixAutomaton<3, 'a'> sa;
        string s;
        for (int i = 0; i < 40; i++){
            s += char('a' + stress::rand_int(0, 1));
            sa.add(s.back());
            int l = stress::rand_int(0, i), r = stress::rand_int(l, i);
            string p = s.substr(l, r - l + 1);
            assert(sa.occurrences(p) == count_occurrences(s, p));
        }
    }

    /// Long inputs: "abbb...b" reaches the 2n - 1 state bound, "aaa...a" is the worst case for occurrences
    const int n = 200000;
    SuffixAutomaton<3, 'a'> skewed('a' + string(n - 1, 'b'));
    assert((int)skewed.len.size() == 2 * n - 1 && skewed.distinct_substrings() == 2 * n - 1);
    SuffixAutomaton<3, 'a'> same(string(n, 'a'));
    assert(same.distinct_substrings() == n && same.occurrences(string(n / 2, 'a')) == n - n / 2 + 1);

    /// Distinct substrings of a long random string against sorted suffixes and their common prefixes
    string big = random_string(3000, 2);
    vector<int> order(big.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y){ return big.compare(x, string::npos, big, y, string::npos) < 0; });
    long long distinct = 0;
    for (int i = 0; i < (int)order.size(); i++){
        int lcp = 0;
        while (i > 0 && max(order[i], order[i - 1]) + lcp < (int)big.size() && big[order[i] + lcp] == big[order[i - 1] + lcp]) lcp++;
        distinct += (long long)big.size() - order[i] - lcp;
    }
    SuffixAutomaton<3, 'a'> sa(big);
    assert(sa.distinct_substrings() == distinct);

    string text = random_string(n, 3);
    SuffixAutomaton<3, 'a'> big_sa(text);
    for (int k = 0; k < 50; k++){
        int l = stress::rand_int(0, n - 1), r = min(n - 1, l + (int)stress::rand_int(0, 8));
        string p = text.substr(l, r - l + 1);
        assert(big_sa.occurrences(p) == count_occurrences(text, p));
    }
    check_lcs(random_string(2000, 2), random_string(2000, 2));

    /// Default 'a'..'z' alphabet at full size, patterns sampled from the text so most counts are nonzero
    string letters = random_string(n, 26);
    SuffixAutomaton<> letters_sa(letters);
    for (int k = 0; k < 50; k++){
        int l = stress::rand_int(0, n - 1), r = min(n - 1, l + (int)stress::rand_int(0, 3));
        string p = letters.substr(l, r - l + 1);
        assert(letters_sa.occurrences(p) == count_occurrences(letters, p));
    }
    assert(letters_sa.occurrences("zz{") == 0 && !letters_sa.contains("A"));

    /// A non-letter alphabet, every check above uses BASE = 'a'
    for (long long it = 0; it < stress::scaled(300); it++){
        string a = random_string(stress::rand_int(0, 12), 2), b = random_string(stress::rand_int(0, 12), 2);
        for (char& ch : a) ch += '0' - 'a';
        for (char& ch : b) ch += '0' - 'a';
        string common = longest_common_substring<2, '0'>(a, b);
        assert((int)common.size() == lcs_length(a, b) && a.find(common) != string::npos && b.find(common) != string::npos);
        SuffixAutomaton<2, '0'> digits(a);
        assert(digits.occurrences(common) == count_occurrences(a, common));
    }

    return 0;
}
