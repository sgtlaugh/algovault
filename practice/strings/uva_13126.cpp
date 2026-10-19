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
    static const uint64_t base;

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

const uint64_t PolyHash::base = mt19937_64(chrono::system_clock::now().time_since_epoch().count())() % (mod / 3) + (mod / 3);

int main(){
    static char text[100005], word[100005];
    while (scanf("%100004s %100004s", text, word) == 2){
        int n = strlen(text), m = strlen(word);
        PolyHash ht(text), hw(word);

        /// At most 10^6 / |t| wildcards, so checking every letter run of w at every start is O(10^6 + |t|)
        vector<pair<int, int>> runs;
        for (int i = 0, j; i < m; i = j){
            for (j = i; j < m && word[j] != '?'; j++) {}
            if (j > i) runs.push_back({i, j - 1});
            j += j == i;
        }
        vector<uint64_t> want;
        for (auto run : runs) want.push_back(hw.get_hash(run.first, run.second));

        int res = 0;
        for (int s = 0; s + m <= n; s++){
            bool ok = true;
            for (int k = 0; k < (int)runs.size() && ok; k++) ok = ht.get_hash(s + runs[k].first, s + runs[k].second) == want[k];
            res += ok;
        }
        printf("%d\n", res);
    }
    return 0;
}
