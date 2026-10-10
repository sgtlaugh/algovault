/***
 *
 * Wildcard Pattern Matching with NTT
 * Finds every position where pattern occurs in text, a wildcard in either string matches any character
 *
 * Complexity: O(n log n) for text length n, three convolutions modulo 998244353
 *
 * wildcard_match(text, pattern, wildcard = '?') returns the 0-indexed start positions in increasing order
 * Works on strings, string literals (mixed freely) and on vectors of any ordered value type, pass the wildcard value for vectors
 * An empty pattern or one longer than the text returns no positions, as kmp_search does
 * Requires n <= 2^23, the largest NTT length modulo 998244353
 * The NTT core is a checked copy of algebra/ntt.cpp
 *
 * Each non-wildcard value gets a random nonzero weight w (0 for the wildcard), and position i scores
 * sum_j w(p[j]) * w(t[i + j]) * (w(p[j]) - w(t[i + j]))^2, which is 0 exactly when every pair matches
 * A real match is always reported, a mismatching position is falsely reported with probability at most 4 / 998244352
 * A whole call reports some false position with probability at most (n - m + 1) * 4 / 998244352 (~3% at n = 2^23)
 * Small alphabets err far less: with two letters every score is a multiple of one term, so a call errs only if both weights collide
 *
 * Example:
 *   wildcard_match("abacaba", "a?a")  // {0, 2, 4}
 *   wildcard_match("a?c?b", "abc")    // {0}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace wildcard_ntt{
    const unsigned MOD = 998244353;
    const unsigned ROOT = 3;
    const int max_log = 23;

    /// BEGIN COPY ntt_core from code_library/algebra/ntt.cpp
    constexpr unsigned long long power(unsigned long long b, unsigned long long e, unsigned long long mod){
        unsigned long long res = 1;
        for (b %= mod; e; e >>= 1, b = b * b % mod){
            if (e & 1) res = res * b % mod;
        }
        return res;
    }

    /// rt[k + j] = w^j for every power of two k < n, where w is a primitive 2k-th root of unity
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n){
        vector<unsigned> rt(max(n, 2), 1);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % MOD : rt[i >> 1];
        }
        return rt;
    }

    /// Forward transform in place, a.size() must be a power of two no larger than rt.size()
    template<unsigned MOD>
    void transform(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++){
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }

        /// MOD < 2^31 keeps u + v below 2^32
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = (unsigned long long)rt[j + k] * a[i + j + k] % MOD;
                    a[i + j] = u + v >= MOD ? u + v - MOD : u + v;
                    a[i + j + k] = u >= v ? u - v : u + MOD - v;
                }
            }
        }
    }
    /// END COPY ntt_core

    /// Transform of w^k placed at the front of a zero-padded array of length len
    vector<unsigned> transformed_powers(const vector<unsigned>& weights, int k, int len, const vector<unsigned>& rt){
        vector<unsigned> a(len);
        for (int i = 0; i < (int)weights.size(); i++) a[i] = power(weights[i], k, MOD);
        transform<MOD>(a, rt);
        return a;
    }
}

template<typename Container>
vector<int> wildcard_match(const Container& text, const Container& pattern, typename Container::value_type wildcard = '?'){
    using namespace wildcard_ntt;
    using T = typename Container::value_type;

    int n = text.size(), m = pattern.size();
    vector<int> positions;
    if (m == 0 || m > n) return positions;
    assert(n <= (1 << max_log));

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    auto random_weight = [&]{ return uniform_int_distribution<unsigned>(1, MOD - 1)(rng); };
    vector<unsigned> wt(n), wp(m);

    /// Equal values share one weight, given out by a byte table or by runs of a sort: a map lookup per element cost 4x the NTTs
    if constexpr (is_integral_v<T> && sizeof(T) == 1){
        array<unsigned, 256> weight_of;
        for (auto& w : weight_of) w = random_weight();
        weight_of[(unsigned char)wildcard] = 0;
        for (int i = 0; i < n; i++) wt[i] = weight_of[(unsigned char)text[i]];
        for (int j = 0; j < m; j++) wp[j] = weight_of[(unsigned char)pattern[m - 1 - j]];
    }
    else{
        vector<pair<T, int>> items;
        items.reserve(n + m);
        for (int i = 0; i < n; i++){
            if (!(text[i] == wildcard)) items.emplace_back(text[i], i);
        }
        for (int j = 0; j < m; j++){
            if (!(pattern[m - 1 - j] == wildcard)) items.emplace_back(pattern[m - 1 - j], n + j);
        }
        sort(items.begin(), items.end(), [](const pair<T, int>& a, const pair<T, int>& b){ return a.first < b.first; });

        unsigned w = 0;
        for (int k = 0; k < (int)items.size(); k++){
            if (k == 0 || items[k - 1].first < items[k].first) w = random_weight();
            int at = items[k].second;
            (at < n ? wt[at] : wp[at - n]) = w;
        }
    }

    /// A cyclic length >= n only wraps the convolution onto indices below m - 1, which are never read
    int len = 1;
    while (len < n) len <<= 1;
    auto rt = roots<MOD, ROOT>(len);
    uint64_t len_inv = power(len, MOD - 2, MOD);

    vector<unsigned> score(len);
    for (int k = 1; k <= 3; k++){
        auto a = transformed_powers(wt, k, len, rt), b = transformed_powers(wp, 4 - k, len, rt);
        uint64_t coef = (k == 2 ? MOD - 2 : 1) * len_inv % MOD;
        for (int i = 0; i < len; i++){
            unsigned& s = score[(len - i) & (len - 1)];
            s = (s + (uint64_t)a[i] * b[i] % MOD * coef) % MOD;
        }
    }
    transform<MOD>(score, rt);

    for (int i = m - 1; i < n; i++){
        if (score[i] == 0) positions.push_back(i - m + 1);
    }

    return positions;
}

vector<int> wildcard_match(const string& text, const string& pattern, char wildcard = '?'){
    return wildcard_match<string>(text, pattern, wildcard);
}

int main(){
    assert((wildcard_match("abacaba", "a?a") == vector<int>{0, 2, 4}));
    assert((wildcard_match("abacaba", "aba") == vector<int>{0, 4}));
    assert((wildcard_match("a?c?b", "abc") == vector<int>{0}));     /// wildcards in the text match too
    assert((wildcard_match("?yz", "*yz", '*') == vector<int>{0}));  /// '?' is a plain letter once '*' is the wildcard
    assert((wildcard_match(string("abacaba"), "c?b") == vector<int>{3}));

    /// Any ordered value type, -1 as the wildcard
    assert((wildcard_match(vector<int>{1, 2, 3, 1, 2}, vector<int>{1, -1}, -1) == vector<int>{0, 3}));
    return 0;
}
