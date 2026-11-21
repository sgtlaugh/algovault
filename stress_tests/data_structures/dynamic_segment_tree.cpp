#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/dynamic_segment_tree.cpp"
#undef main

/// Max prefix and suffix sums under range assignment: a non-commutative combine, fill {0, 0, 0} folds to itself
struct PrefixSums{
    long long sum, pre, suf;
};

struct PrefixAssign{
    static PrefixSums combine(PrefixSums a, PrefixSums b){
        return {a.sum + b.sum, max(a.pre, a.sum + b.pre), max(b.suf, b.sum + a.suf)};
    }
    static PrefixSums apply(PrefixSums, long long x, long long len){
        long long s = x * len;
        return {s, max(0LL, s), max(0LL, s)};
    }
    static long long compose(long long, long long x){ return x; }
};

/// Affine maps x -> a x + b modulo a prime, composed left to right: fill is the identity, which folds to itself
const long long MOD = 998244353;

struct Affine{
    long long a = 1, b = 0;
};

struct ComposeSet{
    static Affine combine(Affine f, Affine g){ return {g.a * f.a % MOD, (g.a * f.b + g.b) % MOD}; }
    static Affine apply(Affine, Affine f, long long){ return f; }
    static Affine compose(Affine, Affine f){ return f; }
};

/// Updates past a node's two halves create at most 4 nodes per level; queries create none
void check_nodes(size_t before, size_t after, long long lo, long long hi){
    assert(after - before <= 4 * (size_t)(__lg(max(1LL, hi - lo)) + 1));
}

/// Small ranges against a plain array, with negative bounds and an add of 0 that must keep a pending lazy
void stress_sum_small(){
    for (long long it = 0; it < stress::scaled(8000); it++){
        long long lo = stress::rand_int(-30, 30), hi = lo + stress::rand_int(0, it % 10 == 0 ? 500 : 40);
        DynamicSegmentTree<long long> st(lo, hi);
        vector<long long> a(hi - lo + 1, 0);

        for (int op = 0; op < 80; op++){
            long long l = stress::rand_int(lo, hi), r = stress::rand_int(l, hi);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(-1000000, 1000000) * stress::rand_int(0, 1);
                size_t before = st.nodes.size();
                st.update(l, r, v);
                check_nodes(before, st.nodes.size(), lo, hi);
                for (long long i = l; i <= r; i++) a[i - lo] += v;
            }
            else{
                size_t before = st.nodes.size();
                long long s = accumulate(a.begin() + (l - lo), a.begin() + (r - lo) + 1, 0LL);
                assert(st.query(l, r) == s);
                assert(st.nodes.size() == before);
            }
        }
    }
}

/// Positions up to +-4e18 against the list of updates, unsigned so lengths near 8e18 wrap exactly as the brute does
void stress_sum_huge(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        long long lo = stress::rand_int(0, 1) ? -4000000000000000000LL : stress::rand_int(-1000, 1000);
        long long hi = stress::rand_int(0, 1) ? 4000000000000000000LL : lo + stress::rand_int(0, 1000000000000LL);
        DynamicSegmentTree<unsigned long long> st(lo, hi);
        vector<array<long long, 3>> updates;

        /// Endpoints near the bounds and near each other, so splits land on both sides of every mid
        auto pick = [&](){
            int side = stress::rand_int(0, 3);
            long long base = side == 0 ? lo : side == 1 ? hi : stress::rand_int(lo, hi);
            long long p = base + stress::rand_int(-3, 3);
            return min(hi, max(lo, p));
        };

        for (int op = 0; op < 60; op++){
            long long l = pick(), r = pick();
            if (l > r) swap(l, r);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(LLONG_MIN, LLONG_MAX);
                size_t before = st.nodes.size();
                st.update(l, r, v);
                check_nodes(before, st.nodes.size(), lo, hi);
                updates.push_back({l, r, v});
            }
            else{
                unsigned long long s = 0;
                for (auto [ul, ur, v] : updates){
                    long long a = max(l, ul), b = min(r, ur);
                    if (a <= b) s += (unsigned long long)v * (unsigned long long)(b - a + 1);
                }
                size_t before = st.nodes.size();
                assert(st.query(l, r) == s);
                assert(st.nodes.size() == before);
            }
        }
    }
}

/// Range add with range min from a random fill: untouched positions and absent children must read as fill
struct MinAdd{
    static long long combine(long long a, long long b){ return min(a, b); }
    static long long apply(long long val, long long lz, long long){ return val + lz; }
    static long long compose(long long a, long long b){ return a + b; }
};

void stress_min_fill(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        long long lo = stress::rand_int(-20, 20), hi = lo + stress::rand_int(0, it % 10 == 0 ? 400 : 30);
        long long fill = stress::rand_int(-1000, 1000);
        DynamicSegmentTree<long long, MinAdd> st(lo, hi, fill);
        vector<long long> a(hi - lo + 1, fill);

        for (int op = 0; op < 60; op++){
            long long l = stress::rand_int(lo, hi), r = stress::rand_int(l, hi);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(-100, 100);
                st.update(l, r, v);
                for (long long i = l; i <= r; i++) a[i - lo] += v;
            }
            else assert(st.query(l, r) == *min_element(a.begin() + (l - lo), a.begin() + (r - lo) + 1));
        }
    }
}

/// Range assignment with a non-commutative combine, which a swapped combine or a lazy applied in the wrong order breaks
void stress_prefix_assign(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        long long lo = stress::rand_int(-20, 20), hi = lo + stress::rand_int(0, it % 10 == 0 ? 400 : 30);
        DynamicSegmentTree<PrefixSums, PrefixAssign, long long> st(lo, hi, PrefixSums{0, 0, 0});
        vector<long long> a(hi - lo + 1, 0);

        for (int op = 0; op < 60; op++){
            long long l = stress::rand_int(lo, hi), r = stress::rand_int(l, hi);
            if (stress::rand_int(0, 1)){
                long long x = stress::rand_int(-20, 20);
                st.update(l, r, x);
                for (long long i = l; i <= r; i++) a[i - lo] = x;
                continue;
            }

            long long s = 0, pre = 0, suf = 0;
            for (long long i = l; i <= r; i++) s += a[i - lo], pre = max(pre, s);
            long long t = 0;
            for (long long i = r; i >= l; i--) t += a[i - lo], suf = max(suf, t);
            PrefixSums v = st.query(l, r);
            assert(v.sum == s && v.pre == pre && v.suf == suf);
        }
    }
}

/// Point set of affine maps folded left to right, the composite judge problem's shape, against a plain array
void stress_compose(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        long long lo = stress::rand_int(0, 5), hi = lo + stress::rand_int(0, it % 10 == 0 ? 300 : 30);
        DynamicSegmentTree<Affine, ComposeSet> st(lo, hi);
        vector<Affine> f(hi - lo + 1);

        for (int op = 0; op < 60; op++){
            if (stress::rand_int(0, 1)){
                long long p = stress::rand_int(lo, hi);
                Affine g{stress::rand_int(1, MOD - 1), stress::rand_int(0, MOD - 1)};
                st.update(p, p, g);
                f[p - lo] = g;
                continue;
            }

            long long l = stress::rand_int(lo, hi), r = stress::rand_int(l, hi);
            Affine acc;
            for (long long i = l; i <= r; i++) acc = ComposeSet::combine(acc, f[i - lo]);
            Affine got = st.query(l, r);
            assert(got.a == acc.a && got.b == acc.b);
        }
    }
}

void check_asserts(){
    stress::expect_assert("data_structures/dynamic_segment_tree.cpp", []{ DynamicSegmentTree<long long> st(5, 4); });
    stress::expect_assert("data_structures/dynamic_segment_tree.cpp", []{
        DynamicSegmentTree<long long> st(0, 10);
        st.query(3, 11);
    });
    stress::expect_assert("data_structures/dynamic_segment_tree.cpp", []{
        DynamicSegmentTree<long long> st(0, 10);
        st.update(-1, 2, 5);
    });
}

int main(){
    stress_sum_small();
    stress_sum_huge();
    stress_min_fill();
    stress_prefix_assign();
    stress_compose();
    check_asserts();

    return 0;
}
