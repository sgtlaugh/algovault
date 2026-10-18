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
    int n;
    scanf("%d", &n);
    vector<int> x(n), y(n), order(n);
    for (int i = 0; i < n; i++) scanf("%d %d", &x[i], &y[i]);

    /// Left end ascending, right end descending: every range that can contain i comes before it
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        return x[a] != x[b] ? x[a] < x[b] : y[a] > y[b];
    });

    vector<int> contains(n), contained(n);
    OrderedMultiset<int> ends;
    for (int i : order){
        contained[i] = ends.size() - ends.count_less(y[i]);
        ends.insert(y[i]);
    }

    ends = OrderedMultiset<int>();
    for (int j = n - 1; j >= 0; j--){
        int i = order[j];
        contains[i] = ends.count_less(y[i] + 1);
        ends.insert(y[i]);
    }

    for (int i = 0; i < n; i++) printf("%d%c", contains[i], i == n - 1 ? '\n' : ' ');
    for (int i = 0; i < n; i++) printf("%d%c", contained[i], i == n - 1 ? '\n' : ' ');
    return 0;
}
