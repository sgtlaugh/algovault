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

int main(){
    const int K = 10;
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        static char buf[1010];
        scanf("%s", buf);
        string s = buf;
        int n = s.size();
        PolyHash h(s);

        /// Ids for the substrings of length at most K, the only ones a query can name
        vector<uint64_t> shorts;
        for (int i = 0; i < n; i++){
            for (int len = 1; len <= K && i + len <= n; len++) shorts.push_back(h.get_hash(i, i + len - 1));
        }
        sort(shorts.begin(), shorts.end());
        shorts.erase(unique(shorts.begin(), shorts.end()), shorts.end());
        auto id_of = [&](uint64_t x){
            auto it = lower_bound(shorts.begin(), shorts.end(), x);
            return it != shorts.end() && *it == x ? int(it - shorts.begin()) : -1;
        };
        vector<array<int, K + 1>> id(n);
        for (int i = 0; i < n; i++){
            for (int len = 1; len <= K && i + len <= n; len++) id[i][len] = id_of(h.get_hash(i, i + len - 1));
        }

        int q;
        scanf("%d", &q);
        vector<int> qx(q), qy(q), res(q, 0);
        vector<vector<int>> asked(shorts.size());
        for (int i = 0; i < q; i++){
            char x[16], y[16];
            scanf("%s %s", x, y);
            qx[i] = id_of(PolyHash(x).get_hash(0, strlen(x) - 1));
            qy[i] = id_of(PolyHash(y).get_hash(0, strlen(y) - 1));
            if (qx[i] != -1 && qy[i] != -1) asked[qx[i]].push_back(i);
        }

        /// Every distinct substring counted once, at its leftmost occurrence, and filed under each of its asked prefixes
        vector<vector<pair<int, int>>> by_prefix(shorts.size());
        vector<pair<uint64_t, int>> windows;
        for (int len = 1; len <= n; len++){
            windows.clear();
            for (int i = 0; i + len <= n; i++) windows.push_back({h.get_hash(i, i + len - 1), i});
            sort(windows.begin(), windows.end());
            for (int w = 0; w < (int)windows.size(); w++){
                if (w && windows[w].first == windows[w - 1].first) continue;
                int i = windows[w].second;
                for (int k = 1; k <= min(K, len); k++){
                    if (!asked[id[i][k]].empty()) by_prefix[id[i][k]].push_back({i, i + len - 1});
                }
            }
        }

        vector<int> cnt(shorts.size(), 0);
        for (int p = 0; p < (int)shorts.size(); p++){
            if (asked[p].empty()) continue;
            for (auto [i, j] : by_prefix[p]){
                for (int l = 1; l <= min(K, j - i + 1); l++) cnt[id[j - l + 1][l]]++;
            }
            for (int qi : asked[p]) res[qi] = cnt[qy[qi]];
            for (auto [i, j] : by_prefix[p]){
                for (int l = 1; l <= min(K, j - i + 1); l++) cnt[id[j - l + 1][l]] = 0;
            }
        }

        printf("Case %d:\n", tc);
        for (int x : res) printf("%d\n", x);
    }
    return 0;
}
