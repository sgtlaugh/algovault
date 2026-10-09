#include "common.h"

#define main library_main
#include "../code_library/lis.cpp"
#undef main

/// res[i] = longest chain ending at i where each earlier element relates to the next by fits(prev, next)
template <typename T, typename F>
vector<int> brute(const vector<T>& a, F fits){
    vector<int> res(a.size(), 1);
    for (size_t i = 0; i < a.size(); i++){
        for (size_t j = 0; j < i; j++) if (fits(a[j], a[i])) res[i] = max(res[i], res[j] + 1);
    }
    return res;
}

int longest(const vector<int>& v){
    return v.empty() ? 0 : *max_element(v.begin(), v.end());
}

template <typename T>
void check(const vector<T>& a){
    auto inc = brute(a, [](const T& x, const T& y){ return x < y; });
    auto non_dec = brute(a, [](const T& x, const T& y){ return x <= y; });
    auto dec = brute(a, [](const T& x, const T& y){ return x > y; });
    auto non_inc = brute(a, [](const T& x, const T& y){ return x >= y; });

    assert(lis_vector(a, false) == inc && lis_vector(a, true) == non_dec);
    assert(lds_vector(a, false) == dec && lds_vector(a, true) == non_inc);
    assert(lis_length(a, false) == longest(inc) && lis_length(a, true) == longest(non_dec));
    assert(lds_length(a, false) == longest(dec) && lds_length(a, true) == longest(non_inc));
    assert(lis_vector(a) == inc && lds_vector(a) == dec && lis_length(a) == longest(inc) && lds_length(a) == longest(dec));  /// defaults are strict
    assert(lis_vector(a, false, greater<T>()) == dec && lis_vector(a, true, greater<T>()) == non_inc);
}

/// A comparator with ties between unequal values, strings ordered by length only
void check_by_length(const vector<string>& s){
    auto shorter = [](const string& x, const string& y){ return x.size() < y.size(); };
    auto not_longer = [](const string& x, const string& y){ return x.size() <= y.size(); };
    assert(lis_vector(s, false, shorter) == brute(s, shorter));
    assert(lis_vector(s, true, shorter) == brute(s, not_longer));
}

int main(){
    for (long long it = 0; it < stress::scaled(1200); it++){
        int n = stress::rand_int(0, it % 10 ? 25 : 400), spread = stress::rand_int(1, it % 3 ? 5 : 1000000000);
        vector<int> a(n);
        vector<string> s(n);
        for (int i = 0; i < n; i++){
            a[i] = stress::rand_int(-spread, spread);
            s[i] = string(stress::rand_int(0, 2), 'a' + stress::rand_int(0, 2));  /// prefixes compare below their extensions
        }
        check(a);
        check(s);
        check_by_length(s);
    }
    return 0;
}
