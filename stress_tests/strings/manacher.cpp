#include "../common.h"

#define main library_main
#include "../../code_library/strings/manacher.cpp"
#undef main

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int alphabet = stress::rand_int(1, it % 5 ? 2 : 26);
        string s(stress::rand_int(0, it % 10 ? 40 : 1000), 'a');
        for (auto& c : s) c = 'a' + stress::rand_int(0, alphabet - 1);

        auto pal = manacher(s);
        int n = s.size();
        assert((int)pal.size() == max(0, 2 * n - 1));
        for (int i = 0; i + 1 < 2 * n; i++){
            /// Center i covers s[l..r], grown outward while the ends match
            int l = i / 2, r = (i + 1) / 2;
            while (l >= 0 && r < n && s[l] == s[r]) l--, r++;
            assert(pal[i] == r - l - 1);
            assert((i - pal[i] + 1) / 2 == l + 1);  /// the header's start formula
        }

        /// The header's recipe: s[l..r] is a palindrome iff pal[l + r] >= r - l + 1
        for (int l = 0; n <= 40 && l < n; l++){
            for (int r = l; r < n; r++){
                bool is_pal = equal(s.begin() + l, s.begin() + r + 1, s.rbegin() + (n - 1 - r));
                assert(is_pal == (pal[l + r] >= r - l + 1));
            }
        }

        vector<long long> v(s.begin(), s.end());
        for (auto& x : v) x = x * 1000000007LL - 5;
        assert(manacher(v) == pal);
    }

    return 0;
}
