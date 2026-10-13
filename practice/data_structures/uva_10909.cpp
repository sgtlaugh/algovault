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

const int MAX = 2000000;

int main(){
    vector<int> survivors;
    for (int x = 1; x <= MAX; x += 2) survivors.push_back(x);

    /// The early rounds strike a large share of the numbers, a linear pass over an array is cheaper for them
    int round = 1;
    for (; survivors[round] < 256; round++){
        int x = survivors[round], k = 0;
        for (int j = 0; j < (int)survivors.size(); j++){
            if ((j + 1) % x) survivors[k++] = survivors[j];
        }
        survivors.resize(k);
    }

    OrderedSet<int> sieve;
    for (int x : survivors) sieve.insert(x);

    /// Strike every x-th survivor from the back, so the positions still to strike do not shift
    for (; round < (int)sieve.size(); round++){
        int x = *sieve.find_by_order(round);
        if (x > (int)sieve.size()) break;
        for (int j = sieve.size() / x * x - 1; j >= 0; j -= x) sieve.erase(sieve.find_by_order(j));
    }

    vector<int> luckies(sieve.begin(), sieve.end());
    vector<bool> lucky(MAX + 1, false);
    for (int x : luckies) lucky[x] = true;

    int n;
    while (scanf("%d", &n) == 1){
        /// Both luckies are odd, and L2 - L1 is minimised by the largest L1 <= n / 2 whose partner is lucky
        int u = -1;
        if (n % 2 == 0){
            for (int i = upper_bound(luckies.begin(), luckies.end(), n / 2) - luckies.begin() - 1; i >= 0; i--){
                if (lucky[n - luckies[i]]){
                    u = luckies[i];
                    break;
                }
            }
        }
        if (u == -1) printf("%d is not the sum of two luckies!\n", n);
        else printf("%d is the sum of %d and %d.\n", n, u, n - u);
    }
    return 0;
}
