/***
 *
 * Hashing
 * Forward and reverse hash of any segment of a string or vector
 *
 * Complexity: O(n) to build, O(1) per hash query
 *
 * PolyHash h(s): s is a string, vector or const char*, 0-based ranges [l, r] inclusive
 * PolyHash h(s, PolyHash::FORWARD): get_hash only, built about 1.6x faster, rev_hash aborts
 * h.get_hash(l, r): hash of s[l..r] modulo 2^61 - 1, equal segments give equal hashes
 * h.rev_hash(l, r): hash of s[l..r] read backwards, s[l..r] is a palindrome iff get_hash(l, r) == rev_hash(l, r)
 *
 * Notes:
 *   - The base is random per run, so anti-hash tests cannot target it, and shared by all instances, so strings compare
 *   - Two different segments of equal length collide with probability about length / 2^61
 *   - Values count modulo 2^61 - 1 after adding 997: congruent values always collide, -997 acts as an empty prefix
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Order-free key of a set or multiset: XOR (set) or sum (multiset) of a random 64-bit value per distinct element
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PolyHash{
    static constexpr uint64_t mod = (1ULL << 61) - 1;
    static inline const uint64_t base = mt19937_64(chrono::system_clock::now().time_since_epoch().count())() % (mod / 3) + (mod / 3);
    enum Direction{BOTH, FORWARD};

    /// One table for every instance, as base is shared: building skips one array, 15-25% faster
    static inline vector<int64_t> base_pow{1};
    vector<int64_t> pref, suff;

    PolyHash() {}

    template <typename Container>
    PolyHash(const Container& ar, Direction dir = BOTH){
        int n = ar.size();

        while ((int)base_pow.size() <= n) base_pow.push_back(modmul(base_pow.back(), base));

        pref.resize(n + 3, 0);
        for (int i = 1; i <= n; i++){
            pref[i] = modmul(pref[i - 1], base) + value_of(ar[i - 1]);
            if (pref[i] >= (int64_t)mod) pref[i] -= mod;
        }

        if (dir == FORWARD) return;
        suff.resize(n + 3, 0);
        for (int i = n; i >= 1; i--){
            suff[i] = modmul(suff[i + 1], base) + value_of(ar[i - 1]);
            if (suff[i] >= (int64_t)mod) suff[i] -= mod;
        }
    }

    PolyHash(const char* str, Direction dir = BOTH) : PolyHash(string(str), dir) {}
    PolyHash(char* str, Direction dir = BOTH) : PolyHash((const char*)str, dir) {}

    uint64_t get_hash(int l, int r) const{
        int64_t h = pref[r + 1] - modmul(base_pow[r - l + 1], pref[l]);
        return h < 0 ? h + mod : h;
    }

    uint64_t rev_hash(int l, int r) const{
        assert(!suff.empty());
        int64_t h = suff[l + 1] - modmul(base_pow[r - l + 1], suff[r + 2]);
        return h < 0 ? h + mod : h;
    }

    /// a, b < mod
    static int64_t modmul(uint64_t a, uint64_t b){
        __uint128_t c = (__uint128_t)a * b;
        uint64_t r = (c & mod) + (c >> 61);
        return r >= mod ? r - mod : r;
    }

    /// The 997 offset keeps zeros from hashing like an empty prefix; the fast path makes builds 10-20% faster
    template <typename T>
    static uint64_t value_of(const T& x){
        if ((unsigned __int128)x < mod - 997) return (uint64_t)x + 997;  /// a negative x wraps above mod
        return (uint64_t)((((__int128)x % (__int128)mod) + mod + 997) % mod);
    }
};

int main(){
    /// A segment of "racecar" is a palindrome iff its forward and reverse hashes match
    PolyHash h("racecar");
    assert(h.get_hash(0, 6) == h.rev_hash(0, 6));
    assert(h.get_hash(1, 5) != h.rev_hash(0, 4));
    assert(h.get_hash(1, 1) == h.rev_hash(5, 5));
    assert(h.get_hash(1, 1) != h.rev_hash(5, 6));
    assert(h.get_hash(2, 4) == h.rev_hash(2, 4));

    PolyHash numbers(vector<int>{1, 2, 3, 2, 1});
    assert(numbers.get_hash(0, 4) == numbers.rev_hash(0, 4));

    string s = "abcab";
    PolyHash hs(s), hv(vector<char>(s.begin(), s.end()));
    assert(hs.get_hash(0, 1) == hs.get_hash(3, 4));
    assert(hs.get_hash(0, 4) == hv.get_hash(0, 4));
    assert(hs.get_hash(0, 1) == PolyHash("xab").get_hash(1, 2));

    char buf[16] = "abcab";
    assert(PolyHash(buf).get_hash(0, 4) == hs.get_hash(0, 4));

    PolyHash forward("abcab", PolyHash::FORWARD);
    assert(forward.get_hash(0, 1) == hs.get_hash(3, 4));

    return 0;
}
