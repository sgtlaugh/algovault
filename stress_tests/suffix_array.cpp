#include "common.h"

#define main library_main
#include "../code_library/suffix_array.cpp"
#undef main

template <typename Container>
void check(const Container& c){
    int n = c.size();
    auto res = suffix_array(c);

    vector<int> expected(n);
    iota(expected.begin(), expected.end(), 0);
    sort(expected.begin(), expected.end(), [&](int a, int b){ return lexicographical_compare(c.begin() + a, c.end(), c.begin() + b, c.end()); });
    assert(res.sa == expected);

    long long distinct = (long long)n * (n + 1) / 2;
    for (int i = 0; i + 1 < n; i++){
        int a = expected[i], b = expected[i + 1], k = 0;
        while (a + k < n && b + k < n && c[a + k] == c[b + k]) k++;
        assert(res.lcp[i] == k);
        distinct -= k;
    }
    assert(res.distinct_substrings() == distinct);
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(0, it % 10 ? 40 : 2000), alphabet = stress::rand_int(1, 4), mode = stress::rand_int(0, 3);
        vector<long long> v(n);
        for (auto& x : v){
            long long s = stress::rand_int(0, alphabet - 1);
            if (mode == 0) x = s;                                  /// zeros, DC3's sentinel value
            else if (mode == 1) x = s - 2;                         /// negatives
            else if (mode == 2) x = s * 1000000000000LL;           /// a range far too large for DC3 buckets
            else x = stress::rand_int(-1000000, 1000000);
        }
        check(v);

        string s(n, 'a');
        for (auto& ch : s) ch = (char)stress::rand_int(0, 2) ? 'a' + stress::rand_int(0, alphabet - 1) : (char)stress::rand_int(-128, 127);
        check(s);
    }
    return 0;
}
