/***
 *
 * Wildcard Pattern Matching with NTT
 * Finds every position where pattern occurs in text, a wildcard in either string matches any character
 *
 * Complexity: O(n log n) for text length n, three convolutions modulo 998244353
 *
 * wildcard_match(text, pattern, wildcard = '?'): the 0-indexed start positions in increasing order
 *
 * Notes:
 *   - Works on strings and string literals mixed freely, and on vectors of any ordered type given the wildcard value
 *   - An empty pattern or one longer than the text returns no positions, as kmp_search does
 *   - Requires n <= 2^23, the largest NTT length modulo 998244353
 *   - A real match is always reported, a mismatching position is reported with probability at most 4 / 998244352
 *   - A call reports some false position with probability at most (n - m + 1) * 4 / 998244352, about 3% at n = 2^23
 *   - Small alphabets err far less: with two letters a call errs only if both random weights collide
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
        b %= mod;
        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }

        return res;
    }

    /// Newton's step y (2 - p y) doubles the correct low bits of 1 / p modulo 2^32, y = p starts with 3 of them
    constexpr unsigned inverse_mod_2_32(unsigned p, unsigned y, int steps){
        return steps ? inverse_mod_2_32(p, y * (2 - p * y), steps - 1) : y;
    }

    /// Values live in [0, BOUND), 2 MOD below 2^30 to skip the final reduction, else MOD, so u - v + BOUND < 2^32
    template<unsigned MOD>
    constexpr unsigned bound(){
        return MOD < (1u << 30) ? 2 * MOD : MOD;
    }

    /// Montgomery reduction, x / 2^32 modulo MOD for x < MOD 2^32, in [0, bound<MOD>()), replacing the slower x % MOD
    template<unsigned MOD>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = 0u - inverse_mod_2_32(MOD, MOD, 4);
        unsigned r = (x + (unsigned long long)((unsigned)x * NEG_INV) * MOD) >> 32;
        return MOD < (1u << 30) || r < MOD ? r : r - MOD;
    }

    /// rt[k + j] = w^j R modulo MOD for each power of two k < n, w a primitive 2k-th root of unity or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), ((unsigned long long)1 << 32) % MOD);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            z = (z << 32) % MOD;
            for (int i = k; i < 2 * k; i++){
                unsigned w = redc<MOD>(rt[i >> 1] * z);
                rt[i] = i & 1 ? (w >= MOD ? w - MOD : w) : rt[i >> 1];
            }
        }

        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two at most rt.size()
    template<unsigned MOD>
    void forward_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = redc<MOD>((unsigned long long)(u - v + BOUND) * rt[j + k]);
                }
            }
        }
    }

    /// Decimation in time with the inverse roots, bit-reversed order in, natural order out, scaled by n
    template<unsigned MOD>
    void inverse_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<MOD>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = u >= v ? u - v : u + BOUND - v;
                }
            }
        }
    }
    /// END COPY ntt_core

    /// Transform of w^k placed at the front of a zero-padded array of length len
    vector<unsigned> transformed_powers(const vector<unsigned>& weights, int k, int len, const vector<unsigned>& rt){
        vector<unsigned> a(len);
        for (int i = 0; i < (int)weights.size(); i++) a[i] = power(weights[i], k, MOD);
        forward_transform<MOD>(a, rt);
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

    /// Equal values share one weight, from a byte table or runs of a sort: a map lookup per element cost 4x the NTTs
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
    auto rt = roots<MOD, ROOT>(len, false);
    uint64_t len_inv = power(len, MOD - 2, MOD);

    vector<unsigned> score(len);
    for (int k = 1; k <= 3; k++){
        auto a = transformed_powers(wt, k, len, rt), b = transformed_powers(wp, 4 - k, len, rt);
        uint64_t coef = (k == 2 ? MOD - 2 : 1) * len_inv % MOD;
        for (int i = 0; i < len; i++){
            score[i] = (score[i] + (uint64_t)a[i] * b[i] % MOD * coef) % MOD;
        }
    }
    inverse_transform<MOD>(score, roots<MOD, ROOT>(len, true));

    for (int i = m - 1; i < n; i++){
        if (score[i] % MOD == 0) positions.push_back(i - m + 1);
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
