#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/treap.cpp"
#undef main

/// Ordered multiset operations against a sorted vector
void check_multiset(int ops, long long range){
    Treap<long long> treap;
    vector<long long> ref;

    for (int op = 0; op < ops; op++){
        long long key = stress::rand_int(-range, range);
        int kind = stress::rand_int(0, 9);
        if (kind < 4){
            treap.insert(key);
            ref.insert(lower_bound(ref.begin(), ref.end(), key), key);
        }
        else if (kind < 7){
            auto it = lower_bound(ref.begin(), ref.end(), key);
            bool present = it != ref.end() && *it == key;
            assert(treap.erase(key) == present);
            if (present) ref.erase(it);
        }

        assert(treap.size() == (int)ref.size());
        assert(treap.count_less(key) == lower_bound(ref.begin(), ref.end(), key) - ref.begin());
        assert(treap.count(key) == upper_bound(ref.begin(), ref.end(), key) - lower_bound(ref.begin(), ref.end(), key));
        if (!ref.empty()){
            int k = stress::rand_int(0, ref.size() - 1);
            assert(treap.kth(k) == ref[k]);
        }
    }
}

/// Sequence operations against a plain vector
void check_sequence(int n, int ops){
    vector<long long> ref(n);
    for (auto& x : ref) x = stress::rand_int(-1000000000, 1000000000);
    ImplicitTreap treap(ref);
    assert(treap.to_vector() == ref);

    for (int op = 0; op < ops; op++){
        int kind = stress::rand_int(0, 8), sz = ref.size();
        if (kind == 0 || sz == 0){
            int pos = stress::rand_int(0, sz);
            long long v = stress::rand_int(-1000000000, 1000000000);
            treap.insert(pos, v);
            ref.insert(ref.begin() + pos, v);
        }
        else if (kind == 1){
            int pos = stress::rand_int(0, sz - 1);
            treap.erase(pos);
            ref.erase(ref.begin() + pos);
        }
        else if (kind == 7){
            /// Empty vectors included: inserting nothing must leave the sequence unchanged
            int pos = stress::rand_int(0, sz);
            vector<long long> values(stress::rand_int(0, 20));
            for (auto& x : values) x = stress::rand_int(-1000000000, 1000000000);
            treap.insert(pos, values);
            ref.insert(ref.begin() + pos, values.begin(), values.end());
        }
        else if (kind == 8){
            int l = stress::rand_int(0, sz - 1), r = stress::rand_int(l, sz - 1);
            treap.erase(l, r);
            ref.erase(ref.begin() + l, ref.begin() + r + 1);
        }
        else{
            int l = stress::rand_int(0, sz - 1), r = stress::rand_int(l, sz - 1);
            if (kind == 2){
                long long delta = stress::rand_int(-1000000, 1000000);
                treap.update(l, r, delta);
                for (int i = l; i <= r; i++) ref[i] += delta;
            }
            else if (kind == 3){
                treap.reverse(l, r);
                std::reverse(ref.begin() + l, ref.begin() + r + 1);
            }
            else if (kind == 4) assert(treap.get(l) == ref[l]);
            else assert(treap.query(l, r) == accumulate(ref.begin() + l, ref.begin() + r + 1, 0LL));
        }

        assert(treap.size() == (int)ref.size());
        if (op % 50 == 0) assert(treap.to_vector() == ref);
    }

    assert(treap.to_vector() == ref);
}

const long long MOD = 1000000007;

/// x -> mul x + add, so compose(old, new) is order sensitive and catches a swapped compose
struct Affine{
    long long mul, add;
};

struct AffineSum{
    static long long combine(long long a, long long b){ return (a + b) % MOD; }
    static long long apply(long long val, Affine lz, int len){ return (val * lz.mul + lz.add * len) % MOD; }
    static Affine compose(Affine old_lz, Affine new_lz){
        return {old_lz.mul * new_lz.mul % MOD, (old_lz.add * new_lz.mul + new_lz.add) % MOD};
    }
    static long long reversed(long long a){ return a; }
};

/// Best subarray sum with range assign: combine is not commutative, so reverse must swap prefix and suffix
struct Seg{
    long long sum, pref, suf, best;
};

struct MaxSubarrayAssign{
    static Seg combine(Seg a, Seg b){
        return {a.sum + b.sum, max(a.pref, a.sum + b.pref), max(b.suf, b.sum + a.suf), max({a.best, b.best, a.suf + b.pref})};
    }
    static Seg apply(Seg, long long x, int len){
        long long all = x * len, top = x > 0 ? all : x;
        return {all, top, top, top};
    }
    static long long compose(long long, long long new_lz){ return new_lz; }
    static Seg reversed(Seg a){ return {a.sum, a.suf, a.pref, a.best}; }
};

Seg leaf(long long x){
    return {x, x, x, x};
}

long long best_subarray(const vector<long long>& a, int l, int r){
    long long res = a[l], cur = 0;
    for (int i = l; i <= r; i++){
        cur = max(cur + a[i], a[i]);
        res = max(res, cur);
    }

    return res;
}

/// Custom Ops through every changing call, against a plain vector
void check_ops(int n, int ops){
    vector<long long> ref(n), sub(n);
    for (auto& x : ref) x = stress::rand_int(0, MOD - 1);
    for (auto& x : sub) x = stress::rand_int(-20, 20);
    ImplicitTreap<long long, AffineSum, Affine> affine(ref);
    vector<Seg> leaves;
    for (long long x : sub) leaves.push_back(leaf(x));
    ImplicitTreap<Seg, MaxSubarrayAssign, long long> best(leaves);

    for (int op = 0; op < ops; op++){
        int sz = ref.size(), kind = stress::rand_int(0, 5);
        if (kind == 0 || sz == 0){
            int pos = stress::rand_int(0, sz);
            long long v = stress::rand_int(0, MOD - 1), w = stress::rand_int(-20, 20);
            affine.insert(pos, v), best.insert(pos, leaf(w));
            ref.insert(ref.begin() + pos, v), sub.insert(sub.begin() + pos, w);
            continue;
        }

        int l = stress::rand_int(0, sz - 1), r = stress::rand_int(l, sz - 1);
        if (kind == 1){
            Affine lz{stress::rand_int(0, MOD - 1), stress::rand_int(0, MOD - 1)};
            affine.update(l, r, lz);
            for (int i = l; i <= r; i++) ref[i] = (ref[i] * lz.mul + lz.add) % MOD;
            long long x = stress::rand_int(-20, 20);
            best.update(l, r, x);
            for (int i = l; i <= r; i++) sub[i] = x;
        }
        else if (kind == 2){
            affine.reverse(l, r), best.reverse(l, r);
            std::reverse(ref.begin() + l, ref.begin() + r + 1), std::reverse(sub.begin() + l, sub.begin() + r + 1);
        }
        else if (kind == 3){
            affine.erase(l, r), best.erase(l, r);
            ref.erase(ref.begin() + l, ref.begin() + r + 1), sub.erase(sub.begin() + l, sub.begin() + r + 1);
        }
        else{
            long long sum = 0;
            for (int i = l; i <= r; i++) sum = (sum + ref[i]) % MOD;
            assert(affine.query(l, r) == sum);
            assert(affine.get(l) == ref[l]);
            assert(best.query(l, r).best == best_subarray(sub, l, r));
        }
    }

    assert(affine.to_vector() == ref);
}

/// max_right on nonnegative sums, where pred(sum <= k) is monotone, against a linear scan
void check_max_right(int n, int ops){
    vector<long long> ref(n);
    for (auto& x : ref) x = stress::rand_int(0, 10);
    ImplicitTreap treap(ref);

    for (int op = 0; op < ops && n > 0; op++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (stress::rand_int(0, 1)){
            long long delta = stress::rand_int(0, 5);
            treap.update(l, r, delta);
            for (int i = l; i <= r; i++) ref[i] += delta;
        }
        else treap.reverse(l, r), std::reverse(ref.begin() + l, ref.begin() + r + 1);

        int from = stress::rand_int(0, n);
        long long k = stress::rand_int(-1, 60), sum = 0;
        int expect = from - 1;
        while (expect + 1 < n && sum + ref[expect + 1] <= k) sum += ref[++expect];
        assert(treap.max_right(from, [&](long long s){ return s <= k; }) == expect);
    }

    assert(treap.to_vector() == ref);
}

template <typename Tree>
int max_depth(const Tree& tree){
    int res = 0;
    vector<pair<int, int>> stack;
    if (tree.root != -1) stack.push_back({tree.root, 1});

    while (!stack.empty()){
        auto [t, d] = stack.back();
        stack.pop_back();
        res = max(res, d);
        for (int c : {tree.nodes[t].l, tree.nodes[t].r}){
            if (c != -1) stack.push_back({c, d + 1});
        }
    }

    return res;
}

/// Keys ordered by the priorities of seed 20260105 build an n deep chain, so a fixed seed fails here
void check_adversarial(int n){
    mt19937 old_rng(20260105);
    vector<unsigned int> priority(n);
    for (auto& p : priority) p = old_rng();

    vector<int> order(n), rank(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int i, int j){ return priority[i] < priority[j]; });
    for (int i = 0; i < n; i++) rank[order[i]] = i;

    Treap<int> set;
    for (int i = 0; i < n; i++) set.insert(rank[i]);
    assert(max_depth(set) <= 200);

    ImplicitTreap seq;
    vector<unsigned int> seen;
    for (int i = 0; i < n; i++){
        int pos = lower_bound(seen.begin(), seen.end(), priority[i]) - seen.begin();
        seen.insert(seen.begin() + pos, priority[i]);
        seq.insert(pos, i);
    }
    assert(max_depth(seq) <= 200);
}

int main(){
    /// kth past the end would read a missing child, and insert past the end would append silently
    stress::expect_assert("data_structures/treap.cpp", []{
        Treap<int> set;
        set.insert(7);
        set.kth(1);
    });
    stress::expect_assert("data_structures/treap.cpp", []{ ImplicitTreap({1, 2}).insert(3, 9); });

    check_adversarial(5000);

    for (long long it = 0; it < stress::scaled(1500); it++){
        check_multiset(stress::rand_int(1, 200), it % 2 ? 5 : 1000000000000000000LL);
        check_sequence(stress::rand_int(0, 40), stress::rand_int(1, 200));
        check_ops(stress::rand_int(0, 40), stress::rand_int(1, 200));
        check_max_right(stress::rand_int(0, 40), stress::rand_int(1, 100));
    }

    /// Large inputs: depth stays logarithmic, so recursion is safe and operations stay fast
    Treap<int> big;
    for (int i = 0; i < 300000; i++) big.insert(i);
    for (int i = 0; i < 300000; i += 1000) assert(big.kth(i) == i && big.count_less(i) == i);

    vector<long long> values(300000);
    iota(values.begin(), values.end(), 0);
    ImplicitTreap seq(values);
    for (int i = 0; i < 2000; i++){
        int l = stress::rand_int(0, 299999), r = stress::rand_int(l, 299999);
        seq.reverse(l, r);
        std::reverse(values.begin() + l, values.begin() + r + 1);
    }
    assert(seq.to_vector() == values);

    return 0;
}
