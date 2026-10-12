#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/coordinate_compression.cpp"
#undef main

template <class T>
void check(const vector<T>& v){
    vector<T> distinct = v, first_seen;
    sort(distinct.begin(), distinct.end());
    distinct.erase(unique(distinct.begin(), distinct.end()), distinct.end());
    for (const T& x : v) if (find(first_seen.begin(), first_seen.end(), x) == first_seen.end()) first_seen.push_back(x);

    vector<T> sorted = v, appearance = v;
    assert(compress(sorted) == distinct);
    assert(compress(appearance, false) == first_seen);
    for (size_t i = 0; i < v.size(); i++){
        assert(sorted[i] == lower_bound(distinct.begin(), distinct.end(), v[i]) - distinct.begin());
        assert(appearance[i] == find(first_seen.begin(), first_seen.end(), v[i]) - first_seen.begin());
    }
}

/// Multiples of the bucket count unordered_map picks after reserve(n) all collide under the identity std::hash
void check_anti_hash(){
    const int n = 100000;
    unordered_map<long long, int> probe;
    probe.reserve(n);
    long long buckets = probe.bucket_count();

    vector<long long> v(n);
    for (int i = 0; i < n; i++) v[i] = buckets * (i / 2 + 1);
    clock_t start = clock();
    auto vals = compress(v, false);
    assert(clock() - start < 2 * CLOCKS_PER_SEC);  /// CPU time, parallel sanitizer jobs inflate wall time
    assert((int)vals.size() == n / 2 && vals[0] == buckets && v[n - 1] == n / 2 - 1);
}

int main(){
    check_anti_hash();

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 30 : 500), spread = stress::rand_int(1, 50);
        vector<int> a(n);
        vector<long long> b(n);
        vector<double> c(n);

        for (int i = 0; i < n; i++){
            a[i] = stress::rand_int(0, 2) ? stress::rand_int(-spread, spread) : (stress::rand_int(0, 1) ? INT_MIN : INT_MAX);
            b[i] = stress::rand_int(-spread, spread) * 1000000000000LL + (stress::rand_int(0, 9) ? 0 : LLONG_MAX / 2);
            c[i] = stress::rand_int(-spread, spread) / 4.0;
            if (c[i] == 0 && stress::rand_int(0, 1)) c[i] = -0.0;  /// -0.0 == 0.0, one id for both
        }

        check(a);
        check(b);
        check(c);
    }

    return 0;
}
