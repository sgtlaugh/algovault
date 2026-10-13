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
    /// Open addressing table of window hashes with their counts, slots from earlier checks are told apart by stamp
    const int SIZE = 1 << 17;
    vector<uint64_t> key(SIZE);
    vector<int> cnt(SIZE), stamp(SIZE, 0);
    int checks = 0;

    int m;
    static char buf[40010];
    while (scanf("%d", &m) == 1 && m){
        scanf("%s", buf);
        string s = buf;
        int n = s.size();
        PolyHash h(s);

        /// Rightmost start of a length len substring seen at least m times, or -1
        auto rightmost = [&](int len){
            checks++;
            int res = -1;
            for (int l = 0; l + len <= n; l++){
                uint64_t x = h.get_hash(l, l + len - 1);
                int k = x & (SIZE - 1);
                while (stamp[k] == checks && key[k] != x) k = (k + 1) & (SIZE - 1);
                if (stamp[k] != checks) stamp[k] = checks, key[k] = x, cnt[k] = 0;
                if (++cnt[k] >= m) res = l;
            }
            return res;
        };

        /// A substring repeated m times keeps its prefixes repeated, so the best length is binary searchable
        int lo = 0, hi = n;
        while (lo < hi){
            int mid = (lo + hi + 1) / 2;
            if (rightmost(mid) != -1) lo = mid;
            else hi = mid - 1;
        }

        if (!lo) puts("none");
        else printf("%d %d\n", lo, rightmost(lo));
    }
    return 0;
}
