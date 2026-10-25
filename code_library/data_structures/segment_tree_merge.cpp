/***
 *
 * Segment Tree Merging
 * Multisets of values in [0, n - 1] with merge, split by key or by rank, and k-th smallest
 *
 * Complexity:
 *   - insert, split_k, split_key, kth, count: O(log n)
 *   - merge: amortized, all merges together cost O(nodes ever created) = O((inserts + splits) log n)
 *
 * SegmentTreeMerge st(n): an empty node pool for sets of values in [0, n - 1], n >= 1
 * st.insert(root, x, c = 1): adds c >= 1 copies of x to the set root, passed by reference and 0 for a new set
 * st.merge(a, b): the union of a and b, which must not be used afterwards, a != b unless both are 0
 * st.split_k(root, k): root keeps its k smallest elements, returns the set of the others, k >= 0
 * st.split_key(root, x): root keeps the elements < x, returns the set of the elements >= x, any int x
 * st.count(root, lo, hi): number of elements in [lo, hi]
 * st.kth(root, k): the k-th smallest element, 0-based, 0 <= k < size(root)
 * st.size(root): number of elements in root
 *
 * Notes:
 *   - A set is a root index, 0 is the empty set, and all sets share one node pool
 *   - Each insert and split creates at most ceil(log2 n) + 1 nodes, and nodes freed by merge are reused first
 *   - To sort ranges of an array (CF 558E), keep one set per sorted run, split at l and r, then merge the runs between
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

    void insert(int& root, int x, long long c = 1){
        assert(0 <= x && x < n && c >= 1);
        root = insert(root, 0, n - 1, x, c);
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

    int merge(int a, int b){
        assert(!a || a != b);
        if (!a || !b) return a ^ b;

        nodes[a].l = merge(nodes[a].l, nodes[b].l);
        nodes[a].r = merge(nodes[a].r, nodes[b].r);
        nodes[a].cnt += nodes[b].cnt;
        free_nodes.push_back(b);
        return a;
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

    long long count(int root, int lo, int hi) const{
        return count(root, 0, n - 1, lo, hi);
    }

    long long count(int cur, int lo, int hi, int l, int r) const{
        if (!cur || r < lo || hi < l) return 0;
        if (l <= lo && hi <= r) return nodes[cur].cnt;

        int mid = lo + (hi - lo) / 2;
        return count(nodes[cur].l, lo, mid, l, r) + count(nodes[cur].r, mid + 1, hi, l, r);
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

    long long size(int root) const{
        return nodes[root].cnt;
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
    /// Two sets merged, then split by rank and by key
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
