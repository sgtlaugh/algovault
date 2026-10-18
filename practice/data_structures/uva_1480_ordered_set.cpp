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

/// Merge sort tree over value ranks storing positions, the k-th smallest in [l, r] descends by counting positions inside
struct MergeSortTree{
    int n;
    vector<vector<int>> tree;
    vector<int> sorted_values;

    MergeSortTree(const vector<int>& values) : n(values.size()), tree(4 * n), sorted_values(values){
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){
            return values[a] < values[b];
        });
        sort(sorted_values.begin(), sorted_values.end());
        build(1, 0, n - 1, order);
    }

    int kth(int l, int r, int k) const{
        int idx = 1, a = 0, b = n - 1;
        while (a != b){
            int m = (a + b) / 2;
            const vector<int>& left = tree[idx << 1];
            int inside = upper_bound(left.begin(), left.end(), r) - lower_bound(left.begin(), left.end(), l);
            if (k < inside) idx = idx << 1, b = m;
            else k -= inside, idx = idx << 1 | 1, a = m + 1;
        }
        return sorted_values[a];
    }

    void build(int idx, int a, int b, const vector<int>& order){
        if (a == b){
            tree[idx] = {order[a]};
            return;
        }
        int m = (a + b) / 2;
        build(idx << 1, a, m, order);
        build(idx << 1 | 1, m + 1, b, order);
        tree[idx].resize(b - a + 1);
        merge(tree[idx << 1].begin(), tree[idx << 1].end(), tree[idx << 1 | 1].begin(), tree[idx << 1 | 1].end(), tree[idx].begin());
    }
};

int main(){
    int ops, case_id = 0;
    char op[16];
    while (scanf("%d", &ops) == 1){
        OrderedSet<int> beads;
        vector<int> chain;
        vector<array<int, 3>> range_queries;
        long long sums[3] = {0, 0, 0};

        while (ops--){
            scanf("%s", op);
            if (op[0] == 'I'){
                int x;
                scanf("%d", &x);
                beads.insert(x);
                chain.push_back(x);
            }
            else if (op[6] == '1'){
                /// Beads only append, so positions 1..t never change and Query_1 can wait for the full chain
                int s, t, k;
                scanf("%d %d %d", &s, &t, &k);
                range_queries.push_back({s - 1, t - 1, k - 1});
            }
            else if (op[6] == '2'){
                int x;
                scanf("%d", &x);
                sums[1] += beads.order_of_key(x) + 1;
            }
            else{
                int k;
                scanf("%d", &k);
                sums[2] += *beads.find_by_order(k - 1);
            }
        }

        if (!range_queries.empty()){
            MergeSortTree mst(chain);
            for (const auto& q : range_queries) sums[0] += mst.kth(q[0], q[1], q[2]);
        }

        printf("Case %d:\n%lld\n%lld\n%lld\n", ++case_id, sums[0], sums[1], sums[2]);
    }
    return 0;
}
