#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

/// Max prefix and suffix sums, empty allowed: a non-commutative merge, so a swapped merge order in descent shows up
struct PrefixSums{
    long long sum, pre, suf;
};

PrefixSums operator+(PrefixSums a, PrefixSums b){
    return {a.sum + b.sum, max(a.pre, a.sum + b.pre), max(b.suf, b.sum + a.suf)};
}

struct Assign{
    long long x;
};

struct AssignBlock{
    long long x;
    int len;
};

/// the default apply is val + lz * len and compose is old + new, these overloads turn both into range assignment
AssignBlock operator*(Assign a, int len){
    return {a.x, len};
}

PrefixSums operator+(PrefixSums, AssignBlock b){
    long long s = b.x * b.len;
    return {s, max(0LL, s), max(0LL, s)};
}

Assign operator+(Assign, Assign b){
    return b;
}

/// The header's arithmetic progression recipe, expressed through the default merge/apply/compose via overloads
struct ApSum{
    long long sum, cnt, idx_sum;
};

struct Linear{
    long long c0, c1;
};

ApSum operator+(ApSum a, ApSum b){
    return {a.sum + b.sum, a.cnt + b.cnt, a.idx_sum + b.idx_sum};
}

Linear operator*(Linear lz, int){
    return lz;
}

ApSum operator+(ApSum v, Linear lz){
    return {v.sum + lz.c0 * v.cnt + lz.c1 * v.idx_sum, v.cnt, v.idx_sum};
}

Linear operator+(Linear a, Linear b){
    return {a.c0 + b.c0, a.c1 + b.c1};
}

int brute_max_right(const vector<long long>& a, int l, long long k){
    long long s = 0;
    int r = l - 1;
    while (r < (int)a.size() && s + a[r] <= k) s += a[r++];
    return r;
}

int brute_min_left(const vector<long long>& a, int r, long long k){
    long long s = 0;
    int l = r + 1;
    while (l > 1 && s + a[l - 2] <= k) s += a[--l - 1];
    return l;
}

int brute_max_right_pre(const vector<long long>& a, int l, long long k){
    long long s = 0, best = 0;
    int r = l - 1;
    while (r < (int)a.size() && max(best, s + a[r]) <= k){
        s += a[r++];
        best = max(best, s);
    }
    return r;
}

int brute_min_left_suf(const vector<long long>& a, int r, long long k){
    long long s = 0, best = 0;
    int l = r + 1;
    while (l > 1 && max(best, s + a[l - 2]) <= k){
        s += a[--l - 1];
        best = max(best, s);
    }
    return l;
}

void stress_default_sum(){
    for (long long it = 0; it < stress::scaled(15000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 2000 : 40);
        vector<long long> a(n, 0);  /// the size constructor starts every element at the identity, 0 for sum
        if (it % 2) for (auto& x : a) x = stress::rand_int(-1000000, 1000000);
        SegmentTree<long long> st = it % 2 ? SegmentTree<long long>(a) : SegmentTree<long long>(n);

        for (int op = 0; op < 100; op++){
            int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(-1000000, 1000000) * stress::rand_int(0, 1);  /// half the updates add 0
                st.update(l, r, v);
                for (int i = l; i <= r; i++) a[i - 1] += v;
            }
            else assert(st.query(l, r) == accumulate(a.begin() + l - 1, a.begin() + r, 0LL));
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
                int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
                long long v = stress::rand_int(0, hi);
                st.update(l, r, v);
                for (int i = l; i <= r; i++) a[i - 1] += v;
                continue;
            }

            long long k = stress::rand_int(0, (long long)hi * min(n, 8) * 2);
            auto pred = [k](long long s){ return s <= k; };
            if (kind == 1){
                int l = stress::rand_int(1, n + 1);
                assert(st.max_right(l, pred) == brute_max_right(a, l, k));
            }
            else{
                int r = stress::rand_int(0, n);
                assert(st.min_left(r, pred) == brute_min_left(a, r, k));
            }
        }
    }
}

/// Descent with a non-commutative merge and range assignment, values of both signs
void stress_prefix_descent(){
    const PrefixSums id = {0, 0, 0};
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 1000 : 30);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-20, 20);
        vector<PrefixSums> init(n);
        for (int i = 0; i < n; i++) init[i] = {a[i], max(0LL, a[i]), max(0LL, a[i])};
        SegmentTree<PrefixSums, Assign> st(init, id);

        for (int op = 0; op < 60; op++){
            int kind = stress::rand_int(0, 2);
            if (kind == 0){
                int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
                long long x = stress::rand_int(-20, 20);
                st.update(l, r, Assign{x});
                for (int i = l; i <= r; i++) a[i - 1] = x;
                continue;
            }

            long long k = stress::rand_int(0, 120);
            if (kind == 1){
                int l = stress::rand_int(1, n + 1);
                assert(st.max_right(l, [k](PrefixSums v){ return v.pre <= k; }) == brute_max_right_pre(a, l, k));
            }
            else{
                int r = stress::rand_int(0, n);
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
            init[i] = {a[i], 1, i + 1};
        }
        SegmentTree<ApSum, Linear> st(init, ApSum{0, 0, 0});

        for (int op = 0; op < 60; op++){
            int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
            if (stress::rand_int(0, 1)){
                long long x = stress::rand_int(-1000, 1000), d = stress::rand_int(-1000, 1000);
                st.update(l, r, Linear{x - d * l, d});
                for (int i = l; i <= r; i++) a[i - 1] += x + d * (i - l);
            }
            else assert(st.query(l, r).sum == accumulate(a.begin() + l - 1, a.begin() + r, 0LL));
        }
    }
}

/// n = 0 through both constructors: build() never reaches a == b on [1, 0], so it must not be called
void check_empty(){
    SegmentTree<long long> from_vec(vector<long long>{}), from_size(0);
    auto any = [](long long){ return true; };
    for (auto* st : {&from_vec, &from_size}){
        assert(st->n == 0);
        assert(st->max_right(1, any) == 0);
        assert(st->min_left(0, any) == 1);
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
