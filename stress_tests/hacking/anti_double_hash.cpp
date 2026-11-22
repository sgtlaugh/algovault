#include "../common.h"

#define main library_main
#include "../../code_library/hacking/anti_double_hash.cpp"
#undef main

long long poly_hash(const string& s, long long base, long long mod){
    long long h = 0;
    for (char c : s) h = ((__int128)h * base + c) % mod;
    return h;
}

/// get_hash(x) minus the real hash of to_string(x) is one constant per mod, x = 1 included (three zero chunks)
void check_offsets(AntiHash& hack, long long base1, long long mod1, long long base2, long long mod2){
    long long offset1 = -1, offset2 = -1;
    vector<unsigned long long> xs = {0, 1, 65535, 1ULL << 48, ~0ULL};
    for (int it = 0; it < 200; it++) xs.push_back(stress::rand_int(0, LLONG_MAX) * 2ULL + stress::rand_int(0, 1));

    for (unsigned long long x : xs){
        unsigned long long h = hack.get_hash(x) ^ 0x5555555555555555ULL;
        string s = hack.to_string(x);
        long long d1 = ((long long)(h >> 32) - poly_hash(s, base1, mod1) + mod1) % mod1;
        long long d2 = ((long long)(h & 0xffffffffULL) - poly_hash(s, base2, mod2) + mod2) % mod2;
        if (offset1 == -1) offset1 = d1, offset2 = d2;
        assert(d1 == offset1 && d2 == offset2);
    }
}

void check(long long base1, long long mod1, long long base2, long long mod2){
    auto hack = make_unique<AntiHash>(base1, mod1, base2, mod2);
    check_offsets(*hack, base1, mod1, base2, mod2);

    string s, t;
    tie(s, t) = hack->solve();
    assert(s.size() == 64 && t.size() == 64 && s != t);
    assert(s.find_first_not_of("01") == string::npos && t.find_first_not_of("01") == string::npos);
    assert(poly_hash(s, base1, mod1) == poly_hash(t, base1, mod1));
    assert(poly_hash(s, base2, mod2) == poly_hash(t, base2, mod2));
}

int main(){
    /// Small mods keep the rho near sqrt(mod1 mod2) steps; mod2 below 1.4e9 keeps solve's start value off every cycle
    for (long long it = 0; it < stress::scaled(4); it++){
        long long mod1 = stress::rand_int(10000, 300000), mod2 = stress::rand_int(10000, 300000);
        check(stress::rand_int(2, mod1 - 1), mod1, stress::rand_int(2, mod2 - 1), mod2);
    }

    /// mod1 above 1431633920 lets the top chunk of a hash be zero, mod1 near 2^32 fills the packed high half
    const long long big[] = {2147483647LL, 4294967291LL, 1431699457LL};
    for (long long mod1 : big){
        long long mod2 = stress::rand_int(50, 200);
        check(stress::rand_int(2, mod1 - 1), mod1, stress::rand_int(2, mod2 - 1), mod2);
    }

    stress::expect_assert("hacking/anti_double_hash.cpp", []{ AntiHash bad(3, 1LL << 32, 5, 1000000007); });
    stress::expect_assert("hacking/anti_double_hash.cpp", []{ AntiHash bad(3, 1000000007, 5, 8589934609LL); });

    return 0;
}
