#include "../common.h"

#define main library_main
#include "../../code_library/hacking/anti_polyonmial_hash.cpp"
#undef main

/// The generated pair differs by a +-1 combination of base powers, so it must collide whatever values the two letters hash to
long long poly_hash(const string& s, long long base, long long mod, long long zero, long long one){
    unsigned __int128 h = 0;
    for (char c : s) h = (h * base + (c == '0' ? zero : one)) % mod;
    return h;
}

int main(){
    srand(stress::seed());
    const long long mods[] = {(1LL << 61) - 1, 1000000007, 998244353, 4294967311LL, 6000000007LL, 1000000000000000003LL};  /// 6000000007: far enough past 2^32 that most residue products overflow 64 bits

    for (long long it = 0; it < stress::scaled(40); it++){
        long long mod = it % 3 ? mods[stress::rand_int(0, 5)] : stress::rand_int(1000000000LL, 1LL << 61);
        long long base = stress::rand_int(2, mod - 1);

        auto ah = make_unique<AntiHash>(base, mod);
        string s1, s2;
        tie(s1, s2) = ah->solve();

        assert(s1.size() == s2.size() && s1 != s2);
        assert(s1.find_first_not_of("01") == string::npos && s2.find_first_not_of("01") == string::npos);
        for (int q = 0; q < 5; q++){
            long long zero = stress::rand_int(0, mod - 1), one = stress::rand_int(0, mod - 1);
            assert(poly_hash(s1, base, mod, zero, one) == poly_hash(s2, base, mod, zero, one));
        }
    }

    return 0;
}
