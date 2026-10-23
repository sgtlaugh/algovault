/***
 *
 * Ordered Set and Multiset
 * Balanced BST that also finds the k-th smallest element and the rank of a value
 *
 * Complexity: O(log n) per operation
 *
 * OrderedSet<T> s: a std::set with two extra calls
 * s.find_by_order(k): iterator to the k-th smallest, 0-based, s.end() if k >= size
 * s.order_of_key(x): how many elements are < x
 * OrderedMultiset<T> m: an empty multiset
 * m.erase(x): removes one copy of x, false if there is none
 * m.insert(x): adds one copy of x
 * m.count(x): how many copies of x
 * m.count_less(x): how many elements are < x
 * m.kth(k): the k-th smallest, 0-based, 0 <= k < size()
 * m.size(): number of elements, copies included
 *
 * Notes:
 *   - Built on the GNU policy-based tree in ext/pb_ds, so it needs libstdc++
 *   - A pb_ds tree with less_equal is not a multiset: equal keys never compare equivalent, so erase and find fail
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

    bool erase(const T& x){
        auto it = t.lower_bound({x, INT_MIN});
        if (it == t.end() || it->first != x) return false;
        t.erase(it);
        return true;
    }

    void insert(const T& x){
        t.insert({x, next_id++});
    }

    int count(const T& x) const{
        return t.order_of_key({x, INT_MAX}) - count_less(x);
    }

    int count_less(const T& x) const{
        return t.order_of_key({x, INT_MIN});
    }

    T kth(int k) const{
        assert(0 <= k && k < size());
        return t.find_by_order(k)->first;
    }

    int size() const{
        return t.size();
    }
};

int main(){
    /// Order statistics on a set
    OrderedSet<int> s;
    for (int x : {5, 1, 9, 5, 3}) s.insert(x);
    assert(s.size() == 4);
    assert(*s.find_by_order(0) == 1);
    assert(*s.find_by_order(3) == 9);
    assert(s.find_by_order(4) == s.end());
    assert(s.order_of_key(5) == 2);
    assert(s.order_of_key(6) == 3);

    /// A multiset keeps the copies of 5
    OrderedMultiset<int> m;
    for (int x : {5, 1, 9, 5, -1, 5}) m.insert(x);
    assert(m.size() == 6);
    assert(m.count(5) == 3);
    assert(m.kth(2) == 5);
    assert(m.kth(5) == 9);
    assert(m.count_less(5) == 2);

    assert(m.erase(5));
    assert(m.count(5) == 2);
    assert(!m.erase(7));
    assert(m.size() == 5);

    return 0;
}
