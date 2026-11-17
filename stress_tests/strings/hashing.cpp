#include "../common.h"

#define main library_main
#include "../../code_library/strings/hashing.cpp"
#undef main

/// Hash equality against segment equality, collisions have probability ~n^2 / 2^61 per array and should never appear
template <typename T>
void check_value_of(T x){
    __int128 r = (__int128)x % (__int128)PolyHash::mod;
    if (r < 0) r += PolyHash::mod;
    assert(PolyHash::value_of(x) == (uint64_t)((r + 997) % PolyHash::mod));
}

int main(){
    stress::expect_assert("strings/hashing.cpp", []{ PolyHash("level", PolyHash::FORWARD).rev_hash(0, 4); });

    /// value_of's fast path ends at mod - 997, where the offset would wrap
    const long long mod = PolyHash::mod;
    for (long long d = -3; d <= 3; d++){
        for (long long x : {mod - 997 + d, mod + d, 2 * mod - 997 + d, -997 + d, d}){
            check_value_of(x);
            check_value_of((unsigned long long)x);
        }
    }
    for (long long x : {LLONG_MIN, LLONG_MAX}) check_value_of(x);
    check_value_of(ULLONG_MAX);
    for (int c = -128; c < 128; c++) check_value_of((char)c), check_value_of((signed char)c);
    for (long long it = 0; it < stress::scaled(100000); it++){
        check_value_of((long long)stress::rng()()), check_value_of((unsigned long long)stress::rng()());
        check_value_of((int)stress::rng()());
    }

    /// modmul must return the canonical residue, a value in [mod, 2 mod) would hash equal segments differently
    for (long long it = 0; it < stress::scaled(1000000); it++){
        const uint64_t mod = PolyHash::mod;
        uint64_t a = stress::rng()() % mod, b = stress::rng()() % mod;
        if (it % 7 == 0) a = mod - 1 - stress::rand_int(0, 3), b = mod - 1 - stress::rand_int(0, 3);
        if (it % 7 == 1) a = stress::rand_int(0, 3), b = it % 2 ? mod - 1 - stress::rand_int(0, 3) : stress::rand_int(0, 3);
        assert((uint64_t)PolyHash::modmul(a, b) == (uint64_t)((unsigned __int128)a * b % mod));
    }

    for (long long it = 0; it < stress::scaled(10000); it++){
        int n = stress::rand_int(1, it % 10 ? 30 : 2000), alphabet = stress::rand_int(1, 3);
        int mode = stress::rand_int(0, 2);
        vector<long long> a(n);
        for (auto& x : a){
            long long v = stress::rand_int(0, alphabet - 1);
            /// below -997, so the 997 offset alone cannot make them nonnegative
            if (mode == 1) x = v - 1000000;
            else if (mode == 2) x = v ? (v == 1 ? LLONG_MIN : LLONG_MAX) : 0;  /// limits of the element type
            else x = v;
        }
        PolyHash h(a), hr(vector<long long>(a.rbegin(), a.rend())), hf(a, PolyHash::FORWARD);
        if (mode == 0){
            string s(a.begin(), a.end());
            PolyHash hs(s);
            /// strings hash like vectors
            assert(hs.get_hash(0, n - 1) == h.get_hash(0, n - 1) && hs.rev_hash(0, n - 1) == h.rev_hash(0, n - 1));
        }

        for (int q = 0; q < 200; q++){
            int len = stress::rand_int(1, n), l1 = stress::rand_int(0, n - len), l2 = stress::rand_int(0, n - len);
            bool equal = std::equal(a.begin() + l1, a.begin() + l1 + len, a.begin() + l2);
            assert((h.get_hash(l1, l1 + len - 1) == h.get_hash(l2, l2 + len - 1)) == equal);
            assert(hf.get_hash(l1, l1 + len - 1) == h.get_hash(l1, l1 + len - 1));
            /// the reversed segment in the reversed array
            assert(h.rev_hash(l1, l1 + len - 1) == hr.get_hash(n - l1 - len, n - l1 - 1));

            /// Different lengths must differ too, without the offset runs of zeros would all hash to 0
            int len2 = stress::rand_int(1, n);
            if (len2 != len) assert(h.get_hash(l1, l1 + len - 1) != h.get_hash(0, len2 - 1));
        }
    }

    /// The shared power table grows under instances built earlier, their hashes must not change
    vector<int> first(5, 7);
    PolyHash early(first);
    uint64_t before = early.get_hash(0, 4);
    PolyHash later(vector<int>(stress::scaled(200000), 7));
    assert(early.get_hash(0, 4) == before && later.get_hash(1000, 1004) == before);

    return 0;
}
