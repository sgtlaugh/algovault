#include "../common.h"

#define main library_main
#include "../../code_library/misc/next_palindrome.cpp"
#undef main

bool is_palindrome(const string& s){
    return equal(s.begin(), s.end(), s.rbegin());
}

int main(){
    /// Counting up from every number below 2e5 until the next palindrome
    for (long long x = 0; x < 200000; x++){
        long long y = x + 1;
        while (!is_palindrome(to_string(y))) y++;
        assert(next_palindrome(to_string(x)) == to_string(y));
    }

    /// Long numbers: the answer must be a palindrome, larger, and no palindrome may lie strictly between
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 60);
        string s(n, '0');
        s[0] = char('1' + stress::rand_int(0, 8));
        for (int i = 1; i < n; i++) s[i] = stress::rand_int(0, 3) ? char('0' + stress::rand_int(0, 9)) : '9';
        string p = next_palindrome(s);
        assert(is_palindrome(p));
        assert(p.size() > s.size() || (p.size() == s.size() && p > s));
        if (p.size() == s.size()){
            /// Any same-length palindrome is fixed by its first half, the one just below p's half must not exceed s
            string half = p.substr(0, (n + 1) / 2);
            int i = half.size() - 1;
            while (i >= 0 && half[i] == '0') half[i--] = '9';
            if (i >= 0 && !(i == 0 && half[0] == '1' && n > 1)){
                half[i]--;
                string below = half + string(half.rbegin() + n % 2, half.rend());
                assert(below <= s);
            }
        }
        else assert(s == string(n, '9'));
    }
    assert(next_palindrome(string(100000, '9')) == "1" + string(99999, '0') + "1");
    return 0;
}
