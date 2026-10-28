#include "../common.h"

#define main library_main
#include "../../code_library/strings/subsequence_automaton.cpp"
#undef main

const long long MODS[] = {1, 2, 7, 998244353, 1LL << 62};

template<typename T>
int brute_next(const vector<T>& s, int i, T c){
    for (int j = i; j < (int)s.size(); j++){
        if (s[j] == c) return j;
    }

    return s.size();
}

template<typename T>
bool brute_is_subsequence(const vector<T>& s, const vector<T>& t){
    size_t k = 0;
    for (size_t j = 0; j < s.size() && k < t.size(); j++) k += s[j] == t[k];

    return k == t.size();
}

/// Suffix recurrence over the distinct next symbols, exact in __int128 for n <= 100, unlike the library's prefix one
template<typename T>
__int128 brute_distinct(const vector<T>& s){
    int n = s.size();
    vector<__int128> after(n + 2, 0);
    for (int i = n - 1; i >= 0; i--){
        set<T> seen(s.begin() + i, s.end());
        for (const T& c : seen) after[i] += 1 + after[brute_next(s, i, c) + 1];
    }

    return after[0];
}

template<typename T>
set<vector<T>> all_subsequences(const vector<T>& s){
    set<vector<T>> res;
    int n = s.size();
    for (int mask = 0; mask < (1 << n); mask++){
        vector<T> t;
        for (int i = 0; i < n; i++){
            if (mask >> i & 1) t.push_back(s[i]);
        }
        res.insert(t);
    }

    return res;
}

/// symbols[0] is never in s, so it covers the absent symbol for both automata
void check_small(const vector<int>& s, const vector<int>& symbols, int alphabet, int first){
    int n = s.size();
    SubsequenceAutomaton dense(s, alphabet, first);
    SparseSubsequenceAutomaton<int> sparse(s);

    for (int i = 0; i <= n + 1; i++){
        for (int c : symbols){
            int expected = brute_next(s, min(i, n), c);
            assert(dense.next(i, c) == expected);
            assert(sparse.next(i, c) == expected);
        }
    }

    auto subsequences = all_subsequences(s);
    assert(distinct_subsequences(s) == (long long)subsequences.size() - 1);
    for (const auto& t : subsequences) assert(dense.is_subsequence(t) && sparse.is_subsequence(t));

    for (int k = 0; k < 30; k++){
        vector<int> t(stress::rand_int(0, n + 2));
        for (int& c : t) c = symbols[stress::rand_int(0, symbols.size() - 1)];
        bool expected = subsequences.count(t);
        assert(dense.is_subsequence(t) == expected);
        assert(sparse.is_subsequence(t) == expected);
    }
}

void check_mods(const vector<long long>& s){
    __int128 exact = brute_distinct(s);
    for (long long mod : MODS) assert(distinct_subsequences(s, mod) == (long long)(exact % mod));
}

int main(){
    /// Every string over {a, b} up to length 10 and over {a, b, c} up to 6, with 'z' as the absent symbol
    for (int sigma = 2; sigma <= 3; sigma++){
        for (int n = 0; n <= (sigma == 2 ? 10 : 6); n++){
            int total = 1;
            for (int i = 0; i < n; i++) total *= sigma;
            for (int code = 0; code < total; code++){
                vector<int> s;
                int x = code;
                for (int i = 0; i < n; i++){
                    s.push_back('a' + x % sigma);
                    x /= sigma;
                }
                vector<int> symbols = {'z'};
                for (int c = 0; c < sigma; c++) symbols.push_back('a' + c);
                check_small(s, symbols, sigma, 'a');
            }
        }
    }

    /// Negative symbols with first < 0, and an alphabet wider than the symbols s uses
    for (long long it = 0; it < stress::scaled(1000); it++){
        int sigma = stress::rand_int(1, 4);
        vector<int> s(stress::rand_int(0, 11));
        for (int& c : s) c = -5 + stress::rand_int(0, sigma - 1);
        vector<int> symbols = {-6};
        for (int c = 0; c < sigma + 1; c++) symbols.push_back(-5 + c);
        check_small(s, symbols, sigma + 1, -5);
    }

    /// Long random strings over alphabets of size 1, 2 and 26 against linear scans
    for (int sigma : {1, 2, 26}){
        for (long long it = 0; it < stress::scaled(20); it++){
            int n = stress::rand_int(0, 3000);
            string s;
            for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, sigma - 1));
            vector<int> v(s.begin(), s.end());
            SubsequenceAutomaton sa(s);

            for (int k = 0; k < 200; k++){
                int i = stress::rand_int(0, n + 1);
                char c = 'a' + stress::rand_int(0, 25);
                assert(sa.next(i, c) == brute_next(v, min(i, n), (int)c));
            }
            for (int k = 0; k < 50; k++){
                string t;
                if (k % 2 == 0 && n > 0){
                    for (int i = 0; i < n; i++){
                        if (stress::rand_int(0, 9) == 0) t += s[i];
                    }
                    if (stress::rand_int(0, 1)) t.insert(t.begin() + stress::rand_int(0, t.size()), 'a' + stress::rand_int(0, sigma));
                }
                else{
                    for (int len = stress::rand_int(0, 60); len > 0; len--) t += char('a' + stress::rand_int(0, sigma - 1));
                }
                assert(sa.is_subsequence(t) == brute_is_subsequence(v, vector<int>(t.begin(), t.end())));
            }
        }
    }

    /// Values up to 1e9 with negatives and LLONG_MIN in the sparse automaton, and distinct counts under every mod
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(0, 100), sigma = stress::rand_int(1, 60);
        vector<long long> pool(sigma);
        for (long long& x : pool) x = stress::rand_int(-1000000000, 1000000000);
        pool[0] = LLONG_MIN;
        vector<long long> s(n);
        for (long long& x : s) x = pool[stress::rand_int(0, sigma - 1)];
        SparseSubsequenceAutomaton<long long> sa(s);

        for (int k = 0; k < 100; k++){
            int i = stress::rand_int(0, n + 1);
            long long c = k % 10 ? pool[stress::rand_int(0, sigma - 1)] : stress::rand_int(-5, 5);
            assert(sa.next(i, c) == brute_next(s, min(i, n), c));

            vector<long long> t(stress::rand_int(0, 4));
            for (long long& x : t) x = pool[stress::rand_int(0, sigma - 1)];
            assert(sa.is_subsequence(t) == brute_is_subsequence(s, t));
        }
        check_mods(s);
    }

    /// 100 distinct symbols give 2^100 - 1, so 1 << 62 catches an overflowing 2 * count
    vector<long long> distinct(100);
    iota(distinct.begin(), distinct.end(), 0);
    check_mods(distinct);

    /// A value 2^32 + 1 must not wrap onto symbol 1
    SubsequenceAutomaton wide(vector<long long>{0, 1}, 2, 0);
    assert(wide.next(0, (1LL << 32) + 1) == 2);
    assert(!wide.is_subsequence(vector<long long>{(1LL << 32) + 1}));

    stress::expect_assert("strings/subsequence_automaton.cpp", []{ SubsequenceAutomaton bad("abz", 3); });

    return 0;
}
