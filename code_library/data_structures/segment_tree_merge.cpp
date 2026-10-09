/***
 *
 * Segment Tree Merging
 * Dynamic segment trees over values [0, n), each one a multiset, with merge and split by key or by rank
 *
 * Complexity: O(log n) per insert, split, kth and count
 *             merge is amortized: all merges together cost O(nodes ever created) = O((inserts + splits) log n)
 *
 * SegmentTreeMerge st(n), values in [0, n), n >= 1
 * A set is a root index, 0 is the empty set, all sets share one node pool
 *   st.insert(root, x, c): adds c >= 1 copies of x to the set root (passed by reference, may start as 0)
 *   st.merge(a, b): returns the union of two different sets (a != b unless both are 0),
 *                   a and b must not be used afterwards, b's shared nodes are recycled
 *   st.split_k(root, k): root keeps its k smallest elements, returns the set of the others
 *   st.split_key(root, x): root keeps the elements < x, returns the set of the elements >= x, any int x
 *   st.kth(root, k): the k-th smallest element, 0-based, needs 0 <= k < size(root)
 *   st.count(root, lo, hi): number of elements in [lo, hi], st.size(root): number of elements
 *
 * Each insert and split creates at most ceil(log2 n) + 1 nodes, nodes freed by merge are reused first
 * Recipe: sorting ranges of an array (CF 558E style) keeps one set per sorted run in a std::set,
 * splits the runs at l and r, then merges the runs in between
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SegmentTreeMerge{
    struct Node{
        int l, r;
        long long cnt;
    };

    int n;
    vector<Node> nodes;  /// nodes[0] is the empty set, its cnt stays 0
    vector<int> free_nodes;

    SegmentTreeMerge(int n) : n(n), nodes(1, Node{0, 0, 0}){
        assert(n >= 1);
    }

    long long count(int root, int lo, int hi) const{
        return count(root, 0, n - 1, lo, hi);
    }

    void insert(int& root, int x, long long c = 1){
        assert(0 <= x && x < n && c >= 1);
        root = insert(root, 0, n - 1, x, c);
    }

    int kth(int root, long long k) const{
        assert(0 <= k && k < size(root));
        int cur = root, lo = 0, hi = n - 1;

        while (lo < hi){
            int mid = lo + (hi - lo) / 2;
            long long left = nodes[nodes[cur].l].cnt;
            if (k < left) cur = nodes[cur].l, hi = mid;
            else k -= left, cur = nodes[cur].r, lo = mid + 1;
        }
        return lo;
    }

    int merge(int a, int b){
        assert(!a || a != b);
        if (!a || !b) return a ^ b;

        nodes[a].l = merge(nodes[a].l, nodes[b].l);
        nodes[a].r = merge(nodes[a].r, nodes[b].r);
        nodes[a].cnt += nodes[b].cnt;
        free_nodes.push_back(b);
        return a;
    }

    long long size(int root) const{
        return nodes[root].cnt;
    }

    int split_k(int& root, long long k){
        assert(k >= 0);
        if (k == 0){
            int rest = root;
            root = 0;
            return rest;
        }
        return split_rest(root, k);
    }

    int split_key(int& root, int x){
        return split_k(root, x > 0 ? count(root, 0, x - 1) : 0);
    }

    long long count(int cur, int lo, int hi, int l, int r) const{
        if (!cur || r < lo || hi < l) return 0;
        if (l <= lo && hi <= r) return nodes[cur].cnt;

        int mid = lo + (hi - lo) / 2;
        return count(nodes[cur].l, lo, mid, l, r) + count(nodes[cur].r, mid + 1, hi, l, r);
    }

    /// child indices are read into locals: new_node may reallocate nodes, so no reference into it survives a call
    int insert(int cur, int lo, int hi, int x, long long c){
        if (!cur) cur = new_node();
        nodes[cur].cnt += c;
        if (lo == hi) return cur;

        int mid = lo + (hi - lo) / 2;
        if (x <= mid){
            int child = insert(nodes[cur].l, lo, mid, x, c);
            nodes[cur].l = child;
        }
        else{
            int child = insert(nodes[cur].r, mid + 1, hi, x, c);
            nodes[cur].r = child;
        }
        return cur;
    }

    int new_node(){
        if (free_nodes.empty()){
            nodes.push_back(Node{0, 0, 0});
            return nodes.size() - 1;
        }

        int cur = free_nodes.back();
        free_nodes.pop_back();
        nodes[cur] = Node{0, 0, 0};
        return cur;
    }

    /// needs k >= 1, so cur never empties and no node is left with cnt 0
    int split_rest(int cur, long long k){
        if (!cur || nodes[cur].cnt <= k) return 0;

        int rest = new_node();
        long long left = nodes[nodes[cur].l].cnt;
        if (k > left){
            int child = split_rest(nodes[cur].r, k - left);
            nodes[rest].r = child;
        }
        else{
            nodes[rest].r = nodes[cur].r;
            nodes[cur].r = 0;
            int child = split_rest(nodes[cur].l, k);
            nodes[rest].l = child;
        }

        nodes[rest].cnt = nodes[cur].cnt - k;
        nodes[cur].cnt = k;
        return rest;
    }
};

int main(){
    SegmentTreeMerge st(10);
    int a = 0, b = 0;
    for (int x : {3, 1, 4, 1, 5}) st.insert(a, x);
    for (int x : {9, 2, 6}) st.insert(b, x);
    assert(st.size(a) == 5);
    assert(st.kth(a, 2) == 3);          /// {1, 1, 3, 4, 5}
    assert(st.count(a, 1, 4) == 4);

    a = st.merge(a, b);                 /// {1, 1, 2, 3, 4, 5, 6, 9}, b is gone
    assert(st.size(a) == 8);
    assert(st.kth(a, 7) == 9);

    int high = st.split_k(a, 3);        /// a = {1, 1, 2}, high = {3, 4, 5, 6, 9}
    assert(st.kth(a, 2) == 2);
    assert(st.kth(high, 0) == 3);

    int top = st.split_key(high, 6);    /// high = {3, 4, 5}, top = {6, 9}
    assert(st.size(high) == 3);
    assert(st.kth(top, 0) == 6);
    return 0;
}
