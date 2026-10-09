#include "../common.h"

#define main library_main
#include "../../code_library/strings/z_algorithm.cpp"
#undef main

template <typename Container>
void check(const Container& c){
    int n = c.size();
    auto z = z_function(c);
    assert((int)z.size() == n);
    for (int i = 0; i < n; i++){
        int k = 0;
        while (i + k < n && c[k] == c[i + k]) k++;
        assert(z[i] == k);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int alphabet = stress::rand_int(1, it % 4 ? 2 : 26);
        string s(stress::rand_int(0, it % 10 ? 40 : 1000), 'a');
        for (auto& c : s) c = 'a' + stress::rand_int(0, alphabet - 1);
        if (it % 5 == 0 && !s.empty()){
            int period = stress::rand_int(1, 5);
            for (size_t i = period; i < s.size(); i++) s[i] = s[i - period];  /// long z boxes that reuse earlier values
        }
        check(s);
        check(vector<int>(s.begin(), s.end()));
    }
    return 0;
}
