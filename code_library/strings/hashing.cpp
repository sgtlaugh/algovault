/***
 *
 * 64-bit hashing for vectors or strings
 * Get the forward and reverse hash of any segment
 * Base is chosen randomly to prevent anti-hash cases from being constructed
 *
 * Complexity - O(n) to build, O(1) for each hash query, O(n) memory per instance
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

constexpr uint64_t mod = (1ULL << 61) - 1;

const uint64_t seed = chrono::system_clock::now().time_since_epoch().count();
const uint64_t base = mt19937_64(seed)() % (mod / 3) + (mod / 3);

/// a, b < mod
int64_t modmul(uint64_t a, uint64_t b){
    __uint128_t c = (__uint128_t)a * b;
    uint64_t r = (c & mod) + (c >> 61);
    return r >= mod ? r - mod : r;
}

struct PolyHash{
    /// Remove suff vector and usage if reverse hash is not required for more speed
    vector<int64_t> pref, suff, base_pow;

    PolyHash() {}

    template <typename T>
    PolyHash(const vector<T>& ar){
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

    /// The 997 offset keeps zeros from hashing like an empty prefix
    template <typename T>
    static uint64_t value_of(const T& x){
        return (uint64_t)((((__int128)x % (__int128)mod) + mod + 997) % mod);
    }

    PolyHash(const char* str)
        : PolyHash(vector<char> (str, str + strlen(str))) {}

    uint64_t get_hash(int l, int r){
        int64_t h = pref[r + 1] - modmul(base_pow[r - l + 1], pref[l]);
        return h < 0 ? h + mod : h;
    }

    uint64_t rev_hash(int l, int r){
        int64_t h = suff[l + 1] - modmul(base_pow[r - l + 1], suff[r + 2]);
        return h < 0 ? h + mod : h;
    }
};

int main(){
    PolyHash H = PolyHash("racecar");

    assert(H.get_hash(0, 6) == H.rev_hash(0, 6));
    assert(H.get_hash(1, 5) != H.rev_hash(0, 4));
    assert(H.get_hash(1, 1) == H.rev_hash(5, 5));
    assert(H.get_hash(1, 1) != H.rev_hash(5, 6));
    assert(H.get_hash(2, 4) == H.rev_hash(2, 4));

    H = PolyHash(vector<int> {1, 2, 3, 2, 1});
    assert(H.get_hash(0, 4) == H.rev_hash(0, 4));

    return 0;
}
