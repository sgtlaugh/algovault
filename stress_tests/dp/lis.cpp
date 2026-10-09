#include "../common.h"

#define main library_main
#include "../../code_library/dp/lis.cpp"
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

/// idx must be a strictly increasing index chain linked by fits, as long as the O(n^2) DP optimum
template <typename T, typename F>
void check_chain(const vector<T>& a, const vector<int>& idx, F fits, const vector<int>& dp){
    assert((int)idx.size() == longest(dp));
    for (size_t k = 0; k < idx.size(); k++){
        assert(0 <= idx[k] && idx[k] < (int)a.size());
        if (k) assert(idx[k - 1] < idx[k] && fits(a[idx[k - 1]], a[idx[k]]));
    }
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

    check_chain(a, lis_indices(a), [](const T& x, const T& y){ return x < y; }, inc);
    check_chain(a, lis_indices(a, true), [](const T& x, const T& y){ return x <= y; }, non_dec);
    check_chain(a, lis_indices(a, false, greater<T>()), [](const T& x, const T& y){ return x > y; }, dec);
    check_chain(a, lis_indices(a, true, greater<T>()), [](const T& x, const T& y){ return x >= y; }, non_inc);
}

/// A comparator with ties between unequal values, strings ordered by length only
void check_by_length(const vector<string>& s){
    auto shorter = [](const string& x, const string& y){ return x.size() < y.size(); };
    auto not_longer = [](const string& x, const string& y){ return x.size() <= y.size(); };

    auto inc = brute(s, shorter), non_dec = brute(s, not_longer);
    assert(lis_vector(s, false, shorter) == inc);
    assert(lis_vector(s, true, shorter) == non_dec);

    check_chain(s, lis_indices(s, false, shorter), shorter, inc);
    check_chain(s, lis_indices(s, true, shorter), not_longer, non_dec);
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

    vector<int> big(200000);
    for (int& x : big) x = stress::rand_int(-1000000000, 1000000000);
    auto idx = lis_indices(big);
    assert((int)idx.size() == lis_length(big));
    for (size_t k = 1; k < idx.size(); k++) assert(idx[k - 1] < idx[k] && big[idx[k - 1]] < big[idx[k]]);

    return 0;
}
