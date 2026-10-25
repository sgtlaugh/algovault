/***
 *
 * Ordered Set and Multiset
 * Balanced BST that also finds the k-th smallest element and the rank of a value
 *
 * Complexity: O(log n) per operation
 *
 * OrderedSet<T> s: a std::set with two extra calls
 * *s.find_by_order(k): the k-th smallest, 0-based, s.end() if k >= size
 * s.order_of_key(x): how many elements are < x
 * OrderedMultiset<T> m: an empty multiset
 * m.insert(x): adds one copy of x
 * m.erase(x): removes one copy of x, false if there is none
 * m.size(): number of elements, copies included
 * m.kth(k): the k-th smallest, 0-based, 0 <= k < size()
 * m.count_less(x): how many elements are < x
 * m.count(x): how many copies of x
 *
 * Built on the GNU policy-based tree in ext/pb_ds, so it needs libstdc++
 * A pb_ds tree with less_equal as the comparator is not a multiset: equal keys never compare equivalent, so erase and find silently fail
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
    int n, k;

    while (scanf("%d %d", &n, &k) == 2 && (n || k)){
        OrderedSet<int> circle;
        for (int i = 1; i <= n; i++) circle.insert(i);

        for (int i = 1, p = 1, len = n; i < n; i++, len--){
            if (i & 1){
                p = (p + k - 1) % len;
                if (p == 0) p = len;
                circle.erase(circle.find_by_order(p - 1));
                if (p > len - 1) p = 1;
            }
            else{
                p = (p - k) % len;
                if (p < 0) p += len;
                circle.erase(circle.find_by_order(p));
                if (p == 0) p = len - 1;
            }
        }

        printf("%d\n", *circle.begin());
    }
    return 0;
}
