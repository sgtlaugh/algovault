#include "../common.h"

#define main library_main
#include "../../code_library/dp/range_lis.cpp"
#undef main

template <typename T>
int brute(const vector<T>& a, int l, int r, bool allow_equal){
    vector<T> tails;
    for (int i = l; i <= r; i++){
        auto it = allow_equal ? upper_bound(tails.begin(), tails.end(), a[i])
                              : lower_bound(tails.begin(), tails.end(), a[i]);
        if (it == tails.end()) tails.push_back(a[i]);
        else *it = a[i];
    }

    return tails.size();
}

/// Every query on a, both modes, empty ranges l = r + 1 included
template <typename T>
void check_all(const vector<T>& a){
    int n = a.size();
    vector<pair<int, int>> queries;
    for (int l = 0; l <= n; l++){
        for (int r = l - 1; r < n; r++) queries.push_back({l, r});
    }
    shuffle(queries.begin(), queries.end(), stress::rng());

    for (bool allow_equal : {false, true}){
        RangeLIS lis(a, allow_equal);
        vector<int> res = lis.query(queries);
        for (size_t i = 0; i < queries.size(); i++){
            assert(res[i] == brute(a, queries[i].first, queries[i].second, allow_equal));
        }
    }
}

template <typename T>
void check_random_queries(const vector<T>& a, int q, bool allow_equal){
    int n = a.size();
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        l = stress::rand_int(0, n - 1), r = stress::rand_int(0, n - 1);
        if (l > r) swap(l, r);
    }

    RangeLIS lis(a, allow_equal);
    vector<int> res = lis.query(queries);
    for (int i = 0; i < q; i++) assert(res[i] == brute(a, queries[i].first, queries[i].second, allow_equal));
}

/// The definition itself: S(c) = min-plus product of S(a) and S(b), read back as a permutation, O(n^3)
vector<int> brute_product(const vector<int>& a, const vector<int>& b){
    int n = a.size();
    auto dist = [&](const vector<int>& p){
        vector<vector<int>> S(n + 1, vector<int>(n + 1, 0));
        for (int i = n - 1; i >= 0; i--) for (int j = 0; j <= n; j++) S[i][j] = S[i + 1][j] + (p[i] < j);
        return S;
    };
    vector<vector<int>> A = dist(a), B = dist(b), C(n + 1, vector<int>(n + 1, INT_MAX));
    for (int i = 0; i <= n; i++) for (int k = 0; k <= n; k++){
        for (int j = 0; j <= n; j++) C[i][k] = min(C[i][k], A[i][j] + B[j][k]);
    }

    vector<int> c(n, -1);
    for (int i = 0; i < n; i++) for (int k = 0; k < n; k++){
        if (C[i][k + 1] - C[i][k] - C[i + 1][k + 1] + C[i + 1][k] == 1) c[i] = k;
    }

    return c;
}

vector<int> permutation(int n){
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    shuffle(p.begin(), p.end(), stress::rng());
    return p;
}

void check_product(){
    UnitMongeProduct product;

    /// Every pair of permutations up to size 5, then random ones up to 30
    for (int n = 0; n <= 5; n++){
        vector<int> a(n);
        iota(a.begin(), a.end(), 0);
        do{
            vector<int> b(n);
            iota(b.begin(), b.end(), 0);
            do assert(product(a, b) == brute_product(a, b)); while (next_permutation(b.begin(), b.end()));
        } while (next_permutation(a.begin(), a.end()));
    }
    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(6, 30);
        vector<int> a = permutation(n), b = permutation(n);
        assert(product(a, b) == brute_product(a, b));
    }

    /// At scale: the identity is neutral and the product is associative
    int n = 3000;
    vector<int> id(n), a = permutation(n), b = permutation(n), c = permutation(n);
    iota(id.begin(), id.end(), 0);
    assert(product(a, id) == a);
    assert(product(id, a) == a);
    assert(product(product(a, b), c) == product(a, product(b, c)));
}

int main(){
    check_product();

    RangeLIS empty(vector<int>{});
    assert((empty.query({{0, -1}}) == vector<int>{0}));

    /// Every permutation up to 7, then every array over {0, 1, 2} up to length 7, where ties decide strict against weak
    for (int n = 1; n <= 7; n++){
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        do check_all(p); while (next_permutation(p.begin(), p.end()));
    }
    for (int n = 1; n <= 7; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 3;
        for (int code = 0; code < total; code++){
            vector<int> a(n);
            for (int i = 0, c = code; i < n; i++, c /= 3) a[i] = c % 3;
            check_all(a);
        }
    }

    /// Random sizes cross the power-of-two splits of the build; values span long long, strings use their own <
    for (int it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(1, 200);
        if (it % 4 == 0) check_all(permutation(min(n, 30)));
        vector<long long> a(n);
        bool few_values = it % 2;
        for (auto& x : a) x = few_values ? stress::rand_int(0, 3) : stress::rand_int(LLONG_MIN, LLONG_MAX);
        check_random_queries(a, 300, it % 4 < 2);
    }
    vector<string> words(150);
    for (auto& w : words) w = string(stress::rand_int(1, 2), 'a' + stress::rand_int(0, 2));
    check_random_queries(words, 500, false);
    check_random_queries(words, 500, true);

    /// Large sorted, reversed and nearly sorted inputs give extreme braids, checked on random queries
    int n = stress::scaled(1) > 1 ? 100000 : 6000;
    vector<int> sorted(n), reversed, nearly;
    iota(sorted.begin(), sorted.end(), 0);
    reversed.assign(sorted.rbegin(), sorted.rend());
    nearly = sorted;
    for (int i = 0; i < n / 100; i++) swap(nearly[stress::rand_int(0, n - 1)], nearly[stress::rand_int(0, n - 1)]);
    for (auto* a : {&sorted, &reversed, &nearly}) check_random_queries(*a, 100, false);
    check_random_queries(permutation(n), 100, false);

    stress::expect_assert("dp/range_lis.cpp", []{ RangeLIS(vector<int>{1, 2, 3}).query({{0, 3}}); });
    stress::expect_assert("dp/range_lis.cpp", []{ RangeLIS(vector<int>{1, 2, 3}).query({{4, 2}}); });

    return 0;
}
