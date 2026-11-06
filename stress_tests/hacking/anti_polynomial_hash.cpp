#include "../common.h"

#define main library_main
#include "../../code_library/hacking/anti_polynomial_hash.cpp"
#undef main

/// The pair differs by a +-1 combination of base powers, so it collides whatever values the two letters hash to
long long poly_hash(const string& s, long long base, long long mod, long long zero, long long one){
    unsigned __int128 h = 0;
    for (char c : s) h = (h * base + (c == '0' ? zero : one)) % mod;
    return h;
}

void check(long long base, long long mod, int n){
    auto ah = make_unique<AntiHash>(base, mod, n);
    assert(ah->len >= 2 * AntiHash::chunk);
    /// The depth the speeds in the header were measured at
    if (n == 2000) assert(ah->lim == 6);
    if (n == 1000) assert(ah->lim == 5);

    string s1, s2;
    tie(s1, s2) = ah->solve();

    assert((int)s1.size() == n && (int)s2.size() == n && s1 != s2);
    assert(s1.find_first_not_of("01") == string::npos && s2.find_first_not_of("01") == string::npos);
    assert(poly_hash(s1, base, mod, 1, 2) == poly_hash(s2, base, mod, 1, 2));
    for (int q = 0; q < 5; q++){
        long long zero = stress::rand_int(0, mod - 1), one = stress::rand_int(0, mod - 1);
        assert(poly_hash(s1, base, mod, zero, one) == poly_hash(s2, base, mod, zero, one));
    }
}

int main(){
    srand(stress::seed());
    /// 6000000007: far enough past 2^32 that most residue products overflow 64 bits
    const long long mods[] = {(1LL << 61) - 1, 1000000007, 998244353, 4294967311LL, 6000000007LL, 1000000000000000003LL};

    for (long long it = 0; it < stress::scaled(40); it++){
        long long mod = it % 3 ? mods[stress::rand_int(0, 5)] : stress::rand_int(1000000000LL, 1LL << 61);
        check(stress::rand_int(2, mod - 1), mod, 2000);
    }

    /// The lengths and moduli of Codeforces Gym 100110F, where b = 2 keeps the low powers unreduced
    const long long gym_mods[] = {99999999999973LL, 100000000000000LL, 99999999999999LL, 1LL << 46};
    for (long long it = 0; it < stress::scaled(60); it++){
        long long mod = it % 2 ? gym_mods[stress::rand_int(0, 3)] : stress::rand_int(4, 100000000000000LL);
        long long base = it % 3 ? stress::rand_int(2, 3) : stress::rand_int(2, mod - 2);
        check(base, mod, it % 4 ? stress::rand_int(1000, 1050) : 1000);
    }

    for (long long it = 0; it < stress::scaled(100); it++){
        int n = stress::rand_int(30, 200);
        long long mod = stress::rand_int(4, n < 60 ? 1000 : 1000000);
        check(stress::rand_int(2, mod - 2), mod, n);
    }

    return 0;
}
