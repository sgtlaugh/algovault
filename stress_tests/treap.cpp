#include "common.h"

#define main library_main
#include "../code_library/treap.cpp"
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
        int kind = stress::rand_int(0, 6), sz = ref.size();
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
        else{
            int l = stress::rand_int(0, sz - 1), r = stress::rand_int(l, sz - 1);
            if (kind == 2){
                long long delta = stress::rand_int(-1000000, 1000000);
                treap.add(l, r, delta);
                for (int i = l; i <= r; i++) ref[i] += delta;
            }
            else if (kind == 3){
                treap.reverse(l, r);
                std::reverse(ref.begin() + l, ref.begin() + r + 1);
            }
            else if (kind == 4) assert(treap.get(l) == ref[l]);
            else assert(treap.sum(l, r) == accumulate(ref.begin() + l, ref.begin() + r + 1, 0LL));
        }
        assert(treap.size() == (int)ref.size());
        if (op % 50 == 0) assert(treap.to_vector() == ref);
    }
    assert(treap.to_vector() == ref);
}

int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        check_multiset(stress::rand_int(1, 200), it % 2 ? 5 : 1000000000000000000LL);
        check_sequence(stress::rand_int(0, 40), stress::rand_int(1, 200));
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
