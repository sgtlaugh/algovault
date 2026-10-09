#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/ordered_set.cpp"
#undef main

/// Multiset and set operations against a sorted vector and std::set
int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long range = it % 2 ? 5 : 1000000000000000000LL;
        OrderedMultiset<long long> m;
        OrderedSet<long long> s;
        vector<long long> ref;
        set<long long> ref_set;
        for (int op = 0; op < 150; op++){
            long long x = stress::rand_int(-range, range);
            if (stress::rand_int(0, 2)){
                m.insert(x), s.insert(x);
                ref.insert(upper_bound(ref.begin(), ref.end(), x), x), ref_set.insert(x);
            }
            else{
                auto pos = lower_bound(ref.begin(), ref.end(), x);
                bool present = pos != ref.end() && *pos == x;
                assert(m.erase(x) == present);
                if (present) ref.erase(pos);
                if (!binary_search(ref.begin(), ref.end(), x)) s.erase(x), ref_set.erase(x);
            }
            assert(m.size() == (int)ref.size() && s.size() == ref_set.size());
            assert(m.count_less(x) == lower_bound(ref.begin(), ref.end(), x) - ref.begin());
            assert(m.count(x) == upper_bound(ref.begin(), ref.end(), x) - lower_bound(ref.begin(), ref.end(), x));
            assert((int)s.order_of_key(x) == distance(ref_set.begin(), ref_set.lower_bound(x)));
            if (!ref.empty()){
                int k = stress::rand_int(0, ref.size() - 1);
                assert(m.kth(k) == ref[k]);
                int j = stress::rand_int(0, ref_set.size() - 1);
                assert(*s.find_by_order(j) == *next(ref_set.begin(), j));
            }
        }
    }
    return 0;
}
