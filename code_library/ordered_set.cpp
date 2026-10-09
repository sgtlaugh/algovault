/***
 *
 * Ordered Set and Multiset (GNU policy based data structures)
 * Balanced BST with order statistics: k-th smallest and rank in O(log n)
 *
 * Complexity: O(log n) per operation
 *
 * OrderedSet<T> s: a std::set with two extras
 *   *s.find_by_order(k): k-th smallest, 0-based (s.end() if k >= size)
 *   s.order_of_key(x): number of elements < x
 *
 * OrderedMultiset<T> m: duplicates allowed
 *   insert(x), erase(x) (removes one copy, false if absent), size()
 *   kth(k): k-th smallest, 0-based, needs 0 <= k < size(), count_less(x), count(x)
 *
 * Do not build the multiset with less_equal as the comparator: erase and find then silently fail
 * because two equal keys never compare equivalent. Pairs of (value, unique id) with less<> avoid that
 *
***/

#include <bits/stdtr1c++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

template <typename T>
using OrderedSet = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;

template <typename T>
struct OrderedMultiset{
    OrderedSet<pair<T, int>> t;
    int next_id = 0;

    void insert(const T& x){
        t.insert({x, next_id++});
    }

    bool erase(const T& x){
        auto it = t.lower_bound({x, INT_MIN});
        if (it == t.end() || it->first != x) return false;
        t.erase(it);
        return true;
    }

    int size() const{
        return t.size();
    }

    T kth(int k) const{
        assert(0 <= k && k < size());
        return t.find_by_order(k)->first;
    }

    int count_less(const T& x) const{
        return t.order_of_key({x, INT_MIN});
    }

    int count(const T& x) const{
        return t.order_of_key({x, INT_MAX}) - count_less(x);
    }
};

int main(){
    OrderedSet<int> s;
    for (int x : {5, 1, 9, 5, 3}) s.insert(x);
    assert(s.size() == 4);
    assert(*s.find_by_order(0) == 1 && *s.find_by_order(3) == 9 && s.find_by_order(4) == s.end());
    assert(s.order_of_key(5) == 2 && s.order_of_key(6) == 3 && s.order_of_key(-1) == 0);

    OrderedMultiset<int> m;
    for (int x : {5, 1, 9, 5, -1, 5}) m.insert(x);
    assert(m.size() == 6 && m.count(5) == 3 && m.count(4) == 0);
    assert(m.kth(0) == -1 && m.kth(2) == 5 && m.kth(4) == 5 && m.kth(5) == 9);
    assert(m.count_less(5) == 2 && m.count_less(10) == 6);
    assert(m.erase(5) && m.count(5) == 2 && m.size() == 5);
    assert(!m.erase(7) && m.size() == 5);
    assert(m.erase(-1) && m.kth(0) == 1);
    return 0;
}
