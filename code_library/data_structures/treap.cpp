/***
 *
 * Treap
 * Randomized balanced binary search tree, as an ordered multiset and as an implicit-key sequence
 *
 * Complexity: expected O(log n) per operation, plus O(k) to build from or insert k values at once
 *
 * Treap<T> set: an empty multiset of keys
 * set.erase(key): removes one copy of key, false if there is none
 * set.insert(key): adds one copy of key
 * set.count(key): how many copies of key
 * set.count_less(key): how many keys are < key
 * set.kth(k): the k-th smallest key, 0-based, 0 <= k < size()
 * set.size(): number of keys, copies included
 *
 * ImplicitTreap<T = long long, Ops = SumAdd<T>, L = T> seq(values): the sequence of values, empty when omitted
 * seq.erase(pos): removes the value at pos, 0 <= pos < size()
 * seq.erase(l, r): removes positions l..r
 * seq.insert(pos, value): inserts value before position pos, 0 <= pos <= size(), pos = size() appends
 * seq.insert(pos, values): inserts the whole vector values before position pos, in order
 * seq.reverse(l, r): reverses the order of positions l..r
 * seq.update(l, r, lz): applies the lazy update lz to positions l..r, adds lz by default
 * seq.get(pos): the value at pos, 0 <= pos < size()
 * seq.max_right(l, pred): largest r with pred(seq.query(l, r)) true, l - 1 if none, 0 <= l <= size()
 * seq.query(l, r): combine over positions l..r, their sum by default
 * seq.size(): length of the sequence
 * seq.to_vector(): the whole sequence in order
 *
 * Notes:
 *   - Positions are 0-based, ranges need 0 <= l <= r < size(), and every aggregate must fit T
 *   - pred must stay false once false as the range grows, it is never called on an empty range
 *   - Ops holds four static functions, like SegmentTree's, L is the lazy type when it differs from T:
 *       - combine(a, b): fold of two adjacent ranges, a on the left
 *       - apply(val, lz, len): value of a range of len elements after lz, len = 1 for a single value
 *       - compose(old, new): one lazy meaning old then new, also called on a stale old whose result is dropped
 *       - reversed(a): the aggregate of the range read backwards, a itself when combine is commutative
 *   - Best subarray sum with reversals: T = {sum, prefix, suffix, best}, reversed swaps prefix and suffix
 *   - Treap is a portable, augmentable alternative to ordered_set.cpp's OrderedMultiset for judges without pb_ds
 *   - Long runs of equal values, up to 4.6e18 positions: run_length_implicit_treap.cpp
 *   - Priorities are seeded from the clock: a fixed seed lets an adversary force an O(n) deep chain
 *   - Erased nodes stay in nodes, so memory grows with every insert: reserve nodes ahead for millions of them
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct Treap{
    struct Node{
        T key;
        unsigned int priority;
        int size, l, r;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    bool erase(const T& key){
        return erase(root, key);
    }

    bool erase(int& t, const T& key){
        if (t == -1) return false;
        if (!(nodes[t].key < key) && !(key < nodes[t].key)){
            t = merge(nodes[t].l, nodes[t].r);
            return true;
        }
        bool found = erase(key < nodes[t].key ? nodes[t].l : nodes[t].r, key);
        if (found) pull(t);
        return found;
    }

    void insert(const T& key){
        nodes.push_back({key, (unsigned int)rng(), 1, -1, -1});
        int a, b;
        split(root, key, a, b);
        root = merge(merge(a, (int)nodes.size() - 1), b);
    }

    int count(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (key < nodes[t].key) t = nodes[t].l;
            else res += size(nodes[t].l) + 1, t = nodes[t].r;
        }

        return res - count_less(key);
    }

    int count_less(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (nodes[t].key < key) res += size(nodes[t].l) + 1, t = nodes[t].r;
            else t = nodes[t].l;
        }

        return res;
    }

    T kth(int k) const{
        assert(0 <= k && k < size());
        for (int t = root; ; ){
            int left = size(nodes[t].l);
            if (k < left) t = nodes[t].l;
            else if (k == left) return nodes[t].key;
            else k -= left + 1, t = nodes[t].r;
        }
    }

    int size() const{
        return size(root);
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            nodes[a].r = merge(nodes[a].r, b);
            pull(a);
            return a;
        }
        nodes[b].l = merge(a, nodes[b].l);
        pull(b);
        return b;
    }

    void pull(int t){
        nodes[t].size = 1 + size(nodes[t].l) + size(nodes[t].r);
    }

    /// Splits t into keys < key and keys >= key
    void split(int t, const T& key, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        if (nodes[t].key < key){
            split(nodes[t].r, key, nodes[t].r, b);
            a = t;
        }
        else{
            split(nodes[t].l, key, a, nodes[t].l);
            b = t;
        }
        pull(t);
    }
};

/// Range add and range sum, the default Ops of ImplicitTreap
template <typename T>
struct SumAdd{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
    static T reversed(T a){ return a; }
};

template <typename T = long long, typename Ops = SumAdd<T>, typename L = T>
struct ImplicitTreap{
    struct Node{
        T value, agg;
        L lazy;
        unsigned int priority;
        int size, l, r;
        bool has_lazy, flip;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    ImplicitTreap(const vector<T>& values = {}){
        root = build(values);
    }

    void erase(int pos){
        erase(pos, pos);
    }

    void erase(int l, int r){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        root = merge(a, c);
    }

    void insert(int pos, const T& value){
        assert(0 <= pos && pos <= size());
        int x = new_node(value);
        root = insert_at(root, pos, x);
    }

    void insert(int pos, const vector<T>& values){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, build(values)), b);
    }

    void reverse(int l, int r){
        on_range(l, r, [&](int t){ apply_flip(t); });
    }

    void update(int l, int r, const L& lz){
        on_range(l, r, [&](int t){ apply_lazy(t, lz); });
    }

    T get(int pos){
        return query(pos, pos);
    }

    template <typename F>
    int max_right(int l, F pred){
        assert(0 <= l && l <= size());
        int a, b, res = l - 1;
        split(root, l, a, b);

        bool empty = true;
        T acc{};
        for (int t = b; t != -1; ){
            push(t);
            int c = nodes[t].l;
            if (c != -1){
                T nxt = empty ? nodes[c].agg : Ops::combine(acc, nodes[c].agg);
                if (!pred(nxt)){
                    t = c;
                    continue;
                }
                acc = nxt, empty = false, res += nodes[c].size;
            }

            T nxt = empty ? nodes[t].value : Ops::combine(acc, nodes[t].value);
            if (!pred(nxt)) break;
            acc = nxt, empty = false, res++;
            t = nodes[t].r;
        }

        root = merge(a, b);
        return res;
    }

    T query(int l, int r){
        assert(0 <= l && l <= r && r < size());
        T res{};
        bool found = false;
        fold(root, 0, l, r, res, found);
        return res;
    }

    int size() const{
        return size(root);
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    vector<T> to_vector(){
        vector<T> res;
        vector<int> stack;

        for (int t = root; t != -1 || !stack.empty(); ){
            for (; t != -1; t = nodes[t].l) push(t), stack.push_back(t);
            t = stack.back();
            stack.pop_back();
            res.push_back(nodes[t].value);
            t = nodes[t].r;
        }

        return res;
    }

    /// Cartesian tree of values on fresh random priorities, built left to right with a stack of the right spine
    int build(const vector<T>& values){
        vector<int> stack;

        for (const T& v : values){
            int t = new_node(v), last = -1;
            while (!stack.empty() && nodes[stack.back()].priority < nodes[t].priority){
                last = stack.back();
                stack.pop_back();
                pull(last);
            }
            nodes[t].l = last;
            if (!stack.empty()) nodes[stack.back()].r = t;
            stack.push_back(t);
        }

        int res = stack.empty() ? -1 : stack[0];
        while (!stack.empty()){
            pull(stack.back());
            stack.pop_back();
        }

        return res;
    }

    int new_node(const T& v){
        nodes.push_back({v, v, L{}, (unsigned int)rng(), 1, -1, -1, false, false});
        return nodes.size() - 1;
    }

    /// Descends until x's priority wins, then splits only that subtree: one pass instead of a split and two merges
    int insert_at(int t, int k, int x){
        if (t == -1) return x;
        if (nodes[x].priority > nodes[t].priority){
            int a, b;
            split(t, k, a, b);
            nodes[x].l = a, nodes[x].r = b;
            pull(x);
            return x;
        }

        push(t);
        int left = size(nodes[t].l);
        if (k <= left){
            int c = insert_at(nodes[t].l, k, x);
            nodes[t].l = c;
        }
        else{
            int c = insert_at(nodes[t].r, k - left - 1, x);
            nodes[t].r = c;
        }
        pull(t);
        return t;
    }

    /// Combines positions l..r of subtree t, whose first position is lo, left to right without restructuring
    void fold(int t, int lo, int l, int r, T& res, bool& found){
        if (t == -1 || r < lo || lo + nodes[t].size - 1 < l) return;
        if (l <= lo && lo + nodes[t].size - 1 <= r){
            res = found ? Ops::combine(res, nodes[t].agg) : nodes[t].agg;
            found = true;
            return;
        }

        push(t);
        int mid = lo + size(nodes[t].l);
        fold(nodes[t].l, lo, l, r, res, found);
        if (l <= mid && mid <= r){
            res = found ? Ops::combine(res, nodes[t].value) : nodes[t].value;
            found = true;
        }
        fold(nodes[t].r, mid + 1, l, r, res, found);
    }

    void pull(int t){
        Node& node = nodes[t];
        T agg = node.value;
        node.size = 1;
        if (node.l != -1) agg = Ops::combine(nodes[node.l].agg, agg), node.size += nodes[node.l].size;
        if (node.r != -1) agg = Ops::combine(agg, nodes[node.r].agg), node.size += nodes[node.r].size;
        node.agg = agg;
    }

    template <typename F>
    void on_range(int l, int r, F f){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        f(b);
        root = merge(merge(a, b), c);
    }

    void apply_lazy(int t, const L& lz){
        if (t == -1) return;
        Node& node = nodes[t];
        node.value = Ops::apply(node.value, lz, 1);
        node.agg = Ops::apply(node.agg, lz, node.size);
        /// compose also runs on a stale lazy and the result is dropped: this select beat a branch by 12%
        L composed = Ops::compose(node.lazy, lz);
        node.lazy = node.has_lazy ? composed : lz;
        node.has_lazy = true;
    }

    /// The children swap at once and agg is read backwards, so a flipped node's own fields are already correct
    void apply_flip(int t){
        if (t == -1) return;
        swap(nodes[t].l, nodes[t].r);
        nodes[t].agg = Ops::reversed(nodes[t].agg);
        nodes[t].flip ^= 1;
    }

    /// Splits t into the first k elements and the rest
    void split(int t, int k, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        push(t);
        if (size(nodes[t].l) < k){
            split(nodes[t].r, k - size(nodes[t].l) - 1, nodes[t].r, b);
            a = t;
        }
        else{
            split(nodes[t].l, k, a, nodes[t].l);
            b = t;
        }
        pull(t);
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            push(a);
            nodes[a].r = merge(nodes[a].r, b);
            pull(a);
            return a;
        }
        push(b);
        nodes[b].l = merge(a, nodes[b].l);
        pull(b);
        return b;
    }

    void push(int t){
        Node& node = nodes[t];
        if (!(node.has_lazy | node.flip)) return;

        int l = node.l, r = node.r;
        if (node.has_lazy){
            apply_lazy(l, node.lazy);
            apply_lazy(r, node.lazy);
            node.has_lazy = false;
        }
        if (node.flip){
            apply_flip(l);
            apply_flip(r);
            node.flip = false;
        }
    }
};

int main(){
    /// Ordered multiset that keeps copies
    Treap<int> set;
    for (int x : {5, 1, 9, 5, -1, 7}) set.insert(x);
    assert(set.kth(2) == 5);                    /// sorted: -1 1 5 5 7 9
    assert(set.count_less(5) == 2);
    assert(set.count(5) == 2);
    assert(set.erase(5));                       /// removes one copy
    assert(set.count(5) == 1);
    assert(!set.erase(4));

    /// Sequence edited by position
    ImplicitTreap seq({1, 2, 3, 4, 5});
    seq.update(1, 3, 10);
    assert((seq.to_vector() == vector<long long>{1, 12, 13, 14, 5}));
    seq.reverse(0, 3);
    assert((seq.to_vector() == vector<long long>{14, 13, 12, 1, 5}));
    seq.insert(2, 100);                         /// before position 2
    seq.erase(0);
    assert((seq.to_vector() == vector<long long>{13, 100, 12, 1, 5}));
    assert(seq.query(1, 3) == 113);
    assert(seq.get(1) == 100);

    /// Whole ranges at once: {7, 8} before position 1, then positions 3..5 (100, 12, 1) removed
    seq.insert(1, vector<long long>{7, 8});
    assert((seq.to_vector() == vector<long long>{13, 7, 8, 100, 12, 1, 5}));
    seq.erase(3, 5);
    assert((seq.to_vector() == vector<long long>{13, 7, 8, 5}));

    /// Last r from position 1 whose running sum stays at most 20: 7 + 8 + 5 = 20, so r = 3
    assert(seq.max_right(1, [](long long s){ return s <= 20; }) == 3);

    /// Range assign and range min through an Ops struct, min reads the same both ways
    struct MinAssign{
        static int combine(int a, int b){ return min(a, b); }
        static int apply(int, int lz, int){ return lz; }
        static int compose(int, int lz){ return lz; }
        static int reversed(int a){ return a; }
    };

    ImplicitTreap<int, MinAssign> mins({4, 2, 6, 3});
    mins.update(0, 1, 9);                       /// 9 9 6 3
    assert(mins.query(0, 2) == 6);
    mins.reverse(1, 3);                         /// 9 3 6 9
    assert(mins.get(1) == 3);

    return 0;
}
