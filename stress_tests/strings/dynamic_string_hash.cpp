#include "../common.h"

#define main library_main
#include "../../code_library/strings/dynamic_string_hash.cpp"
#undef main

/// Direct polynomial hash with the same base, sum of (c + 1) * base^(len - 1 - i)
unsigned long long direct(const string& s, int l, int r, unsigned long long base){
    unsigned __int128 h = 0;
    for (int i = l; i <= r; i++) h = (h * base + (unsigned char)s[i] + 1) % DynamicStringHash::MOD;
    return (unsigned long long)h;
}

void check(int n, int ops, int sigma){
    string s(n, 'a');
    for (auto& c : s) c = stress::rand_int(0, 19) ? char('a' + stress::rand_int(0, sigma - 1)) : char(stress::rand_int(-128, 127));
    DynamicStringHash h(s, stress::rng()());

    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (stress::rand_int(0, 2) == 0){
            char c = stress::rand_int(0, 9) ? char('a' + stress::rand_int(0, sigma - 1)) : char(stress::rand_int(-128, 127));
            h.assign(l, r, c);
            fill(s.begin() + l, s.begin() + r + 1, c);
        }
        assert(h.hash(l, r) == direct(s, l, r, h.base));

        int len = r - l + 1, l2 = stress::rand_int(0, n - len);
        bool same = s.compare(l, len, s, l2, len) == 0;
        assert((h.hash(l, r) == h.hash(l2, l2 + len - 1)) == same);
    }
}

/// Instances built with the default seed must share a base, otherwise their hashes cannot be compared
void check_default_seed(){
    DynamicStringHash first("margherita"), second("pepperoni margherita");
    assert(first.base == second.base);
    assert(first.hash(0, 9) == second.hash(10, 19));
    second.assign(0, 8, 'z');
    assert(first.hash(0, 9) == second.hash(10, 19));
}

int main(){
    check_default_seed();
    for (long long it = 0; it < stress::scaled(3000); it++) check(stress::rand_int(1, 50), 100, stress::rand_int(1, 3));
    for (int n : {1, 2, 3, 63, 64, 65}) check(n, 3000, 2);
    check(200000, 2000, 26);
    return 0;
}
