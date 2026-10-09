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

char random_char(int sigma, int odd_one_in){
    return stress::rand_int(0, odd_one_in - 1) ? char('a' + stress::rand_int(0, sigma - 1)) : char(stress::rand_int(-128, 127));
}

/// Short ranges half the time so palindromes show up even in long random strings
pair<int, int> random_range(int n){
    int len = stress::rand_int(0, 1) ? stress::rand_int(1, min(n, 6)) : stress::rand_int(1, n);
    int l = stress::rand_int(0, n - len);
    return {l, l + len - 1};
}

/// Hashes against direct evaluation, equality and palindrome answers against direct string comparison
void check(int n, int ops, int sigma){
    string s(n, 'a');
    for (auto& c : s) c = random_char(sigma, 20);
    DynamicStringHash h(s, stress::rng()());

    for (int op = 0; op < ops; op++){
        auto [l, r] = random_range(n);
        if (stress::rand_int(0, 2) == 0){
            char c = random_char(sigma, 10);
            if (stress::rand_int(0, 1)) r = l;
            h.assign(l, r, c);
            fill(s.begin() + l, s.begin() + r + 1, c);
        }

        string sub = s.substr(l, r - l + 1), reversed(sub.rbegin(), sub.rend());
        assert(h.hash(l, r) == direct(s, l, r, h.base));
        assert(h.rev_hash(l, r) == direct(reversed, 0, r - l, h.base));
        assert(h.is_palindrome(l, r) == (sub == reversed));

        int len = r - l + 1, l2 = stress::rand_int(0, n - len);
        bool same = s.compare(l2, len, sub) == 0, same_reversed = s.compare(l2, len, reversed) == 0;
        assert((h.hash(l, r) == h.hash(l2, l2 + len - 1)) == same);
        assert((h.rev_hash(l, r) == h.hash(l2, l2 + len - 1)) == same_reversed);
    }
}

/// Every range of every binary string up to length 10, after one random assignment
void check_exhaustive(){
    for (int n = 1; n <= 10; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s(n, 'a');
            for (int i = 0; i < n; i++) s[i] = 'a' + (mask >> i & 1);
            DynamicStringHash h(s, stress::rng()());
            auto [l, r] = random_range(n);
            char c = 'a' + stress::rand_int(0, 1);
            h.assign(l, r, c);
            fill(s.begin() + l, s.begin() + r + 1, c);

            for (int i = 0; i < n; i++){
                for (int j = i; j < n; j++){
                    string sub = s.substr(i, j - i + 1);
                    assert(h.is_palindrome(i, j) == equal(sub.begin(), sub.end(), sub.rbegin()));
                    assert(h.rev_hash(i, j) == direct(string(sub.rbegin(), sub.rend()), 0, j - i, h.base));
                }
            }
        }
    }
}

/// Instances built with the default seed must share a base, otherwise their hashes cannot be compared
void check_default_seed(){
    DynamicStringHash first("margherita"), second("pepperoni margherita"), third("atirehgram");
    assert(first.base == second.base && first.base == third.base);
    assert(first.hash(0, 9) == second.hash(10, 19));
    assert(first.hash(0, 9) == third.rev_hash(0, 9));
    second.assign(0, 8, 'z');
    assert(first.hash(0, 9) == second.hash(10, 19));
}

int main(){
    check_default_seed();
    check_exhaustive();

    for (long long it = 0; it < stress::scaled(2000); it++) check(stress::rand_int(1, 50), 100, stress::rand_int(1, 3));
    for (int n : {1, 2, 3, 63, 64, 65}) check(n, 3000, 2);
    check(200000, 2000, 26);
    check(200000, 500, 1);

    return 0;
}
