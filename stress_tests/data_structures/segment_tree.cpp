#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

/// Max prefix and suffix sums, empty allowed: a non-commutative combine, so a swapped combine order in descent shows up
struct PrefixSums{
    long long sum, pre, suf;
};

struct PrefixAssign{
    static PrefixSums combine(PrefixSums a, PrefixSums b){
        return {a.sum + b.sum, max(a.pre, a.sum + b.pre), max(b.suf, b.sum + a.suf)};
    }
    static PrefixSums apply(PrefixSums, long long x, int len){
        long long s = x * len;
        return {s, max(0LL, s), max(0LL, s)};
    }
    static long long compose(long long, long long x){ return x; }
};

/// The header's arithmetic progression recipe
struct ApSum{
    long long sum, cnt, idx_sum;
};

struct Linear{
    long long c0, c1;
};

struct ProgressionAdd{
    static ApSum combine(ApSum a, ApSum b){ return {a.sum + b.sum, a.cnt + b.cnt, a.idx_sum + b.idx_sum}; }
    static ApSum apply(ApSum v, Linear lz, int){ return {v.sum + lz.c0 * v.cnt + lz.c1 * v.idx_sum, v.cnt, v.idx_sum}; }
    static Linear compose(Linear a, Linear b){ return {a.c0 + b.c0, a.c1 + b.c1}; }
};

int brute_max_right(const vector<long long>& a, int l, long long k){
    long long s = 0;
    int r = l - 1;
    while (r + 1 < (int)a.size() && s + a[r + 1] <= k) s += a[++r];
    return r;
}

int brute_min_left(const vector<long long>& a, int r, long long k){
    long long s = 0;
    int l = r + 1;
    while (l > 0 && s + a[l - 1] <= k) s += a[--l];
    return l;
}

int brute_max_right_pre(const vector<long long>& a, int l, long long k){
    long long s = 0, best = 0;
    int r = l - 1;
    while (r + 1 < (int)a.size() && max(best, s + a[r + 1]) <= k){
        s += a[++r];
        best = max(best, s);
    }
    return r;
}

int brute_min_left_suf(const vector<long long>& a, int r, long long k){
    long long s = 0, best = 0;
    int l = r + 1;
    while (l > 0 && max(best, s + a[l - 1]) <= k){
        s += a[--l];
        best = max(best, s);
    }
    return l;
}

void stress_default_sum(){
    for (long long it = 0; it < stress::scaled(15000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 2000 : 40);
        long long fill = stress::rand_int(-1000000, 1000000) * stress::rand_int(0, 1);
        vector<long long> a(n, fill);   /// the size constructor builds n copies of its value, not only of 0
        if (it % 2) for (auto& x : a) x = stress::rand_int(-1000000, 1000000);
        SegmentTree<long long> st = it % 2 ? SegmentTree<long long>(a) : SegmentTree<long long>(n, fill);

        for (int op = 0; op < 100; op++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(-1000000, 1000000) * stress::rand_int(0, 1);  /// half the updates add 0
                st.update(l, r, v);
                for (int i = l; i <= r; i++) a[i] += v;
            }
            else assert(st.query(l, r) == accumulate(a.begin() + l, a.begin() + r + 1, 0LL));
        }
    }
}

/// max_right and min_left on sums of non-negative values (a monotone predicate) against linear prefix scans
void stress_sum_descent(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 1000 : 30);
        int hi = stress::rand_int(0, 1) ? 3 : 1000;
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(0, hi);
        SegmentTree<long long> st(a);

        for (int op = 0; op < 60; op++){
            int kind = stress::rand_int(0, 2);
            if (kind == 0){
                int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
                long long v = stress::rand_int(0, hi);
                st.update(l, r, v);
                for (int i = l; i <= r; i++) a[i] += v;
                continue;
            }

            long long k = stress::rand_int(0, (long long)hi * min(n, 8) * 2);
            auto pred = [k](long long s){ return s <= k; };
            if (kind == 1){
                int l = stress::rand_int(0, n);
                assert(st.max_right(l, pred) == brute_max_right(a, l, k));
            }
            else{
                int r = stress::rand_int(-1, n - 1);
                assert(st.min_left(r, pred) == brute_min_left(a, r, k));
            }
        }
    }
}

/// Descent with a non-commutative combine and range assignment, values of both signs
void stress_prefix_descent(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 1000 : 30);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-20, 20);
        vector<PrefixSums> init(n);
        for (int i = 0; i < n; i++) init[i] = {a[i], max(0LL, a[i]), max(0LL, a[i])};
        SegmentTree<PrefixSums, PrefixAssign, long long> st(init);    /// no identity passed, descent must not need one

        for (int op = 0; op < 60; op++){
            int kind = stress::rand_int(0, 2);
            if (kind == 0){
                int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
                long long x = stress::rand_int(-20, 20);
                st.update(l, r, x);
                for (int i = l; i <= r; i++) a[i] = x;
                continue;
            }

            long long k = stress::rand_int(0, 120);
            if (kind == 1){
                int l = stress::rand_int(0, n);
                assert(st.max_right(l, [k](PrefixSums v){ return v.pre <= k; }) == brute_max_right_pre(a, l, k));
            }
            else{
                int r = stress::rand_int(-1, n - 1);
                assert(st.min_left(r, [k](PrefixSums v){ return v.suf <= k; }) == brute_min_left_suf(a, r, k));
            }
        }
    }
}

/// The arithmetic progression payload from the header against a plain array
void stress_progression_add(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 1000 : 30);
        vector<long long> a(n);
        vector<ApSum> init(n);
        for (int i = 0; i < n; i++){
            a[i] = stress::rand_int(-1000, 1000);
            init[i] = {a[i], 1, i};
        }
        SegmentTree<ApSum, ProgressionAdd, Linear> st(init, ApSum{0, 0, 0});

        for (int op = 0; op < 60; op++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            if (stress::rand_int(0, 1)){
                long long x = stress::rand_int(-1000, 1000), d = stress::rand_int(-1000, 1000);
                st.update(l, r, Linear{x - d * l, d});
                for (int i = l; i <= r; i++) a[i] += x + d * (i - l);
            }
            else assert(st.query(l, r).sum == accumulate(a.begin() + l, a.begin() + r + 1, 0LL));
        }
    }
}

/// n = 0 through both constructors: build() never reaches a == b on [0, -1], so it must not be called
void check_empty(){
    SegmentTree<long long> from_vec(vector<long long>{}), from_size(0);
    auto any = [](long long){ return true; };
    for (auto* st : {&from_vec, &from_size}){
        assert(st->n == 0);
        assert(st->max_right(0, any) == -1);
        assert(st->min_left(-1, any) == 0);
    }
}

int main(){
    check_empty();
    stress_default_sum();
    stress_sum_descent();
    stress_prefix_descent();
    stress_progression_add();

    return 0;
}
