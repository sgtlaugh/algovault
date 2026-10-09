#include "../common.h"

#define main library_main
#include "../../code_library/strings/minimum_rotation.cpp"
#undef main

/// Smallest index of a least rotation by comparing all n rotations
template <typename Container>
int brute(const Container& s){
    int n = s.size(), best = 0;
    auto rotation = [&](int i){
        Container r(s.begin() + i, s.end());
        r.insert(r.end(), s.begin(), s.begin() + i);
        return r;
    };

    for (int i = 1; i < n; i++){
        if (rotation(i) < rotation(best)) best = i;
    }

    return best;
}

int main(){
    for (long long it = 0; it < stress::scaled(30000); it++){
        int n = stress::rand_int(0, 30), sigma = stress::rand_int(1, 3);
        string s;
        for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, sigma - 1));
        if (it % 5 == 0 && n){
            string block = s.substr(0, stress::rand_int(1, n));
            s.clear();
            while ((int)s.size() < n) s += block;
        }
        assert(minimum_rotation(s) == brute(s));

        vector<int> v(n);
        for (auto& x : v) x = stress::rand_int(-2, 2);
        assert(minimum_rotation(v) == brute(v));
    }

    for (int n = 0; n <= 14; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s;
            for (int i = 0; i < n; i++) s += mask >> i & 1 ? 'b' : 'a';
            assert(minimum_rotation(s) == brute(s));
        }
    }

    string big(500000, 'a');
    big[123456] = 'b';
    assert(minimum_rotation(big) == 123457);

    return 0;
}
