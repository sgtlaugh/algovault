/***
 *
 * Treap
 * Randomized balanced binary search tree, as an ordered multiset and as an implicit-key sequence
 *
 * Complexity: expected O(log n) per operation, O(n) to build an ImplicitTreap from an array
 *
 * Treap<T>: ordered multiset of keys
 *   insert(key), erase(key) (removes one copy, false if absent), size()
 *   count_less(key): number of keys < key, count(key): copies of key
 *   kth(k): the k-th smallest key, 0-based, needs 0 <= k < size()
 *
 * ImplicitTreap: sequence of long long values, 0-based positions, ranges [l, r] inclusive
 *   ImplicitTreap(values), insert(pos, value) (before position pos, pos = size() appends), erase(pos), size()
 *   add(l, r, delta), sum(l, r), reverse(l, r), get(pos), to_vector()
 *   Sums must fit in long long
 *
 * Treap is a portable, augmentable alternative to ordered_set.cpp's OrderedMultiset for judges without pb_ds
 * Priorities are seeded from the clock, a fixed seed lets an adversary force an O(n) deep chain
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

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int size() const{
        return size(root);
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

    void insert(const T& key){
        nodes.push_back({key, (unsigned int)rng(), 1, -1, -1});
        int a, b;
        split(root, key, a, b);
        root = merge(merge(a, (int)nodes.size() - 1), b);
    }

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

    int count_less(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (nodes[t].key < key) res += size(nodes[t].l) + 1, t = nodes[t].r;
            else t = nodes[t].l;
        }
        return res;
    }

    int count(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (key < nodes[t].key) t = nodes[t].l;
            else res += size(nodes[t].l) + 1, t = nodes[t].r;
        }
        return res - count_less(key);
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
};

struct ImplicitTreap{
    struct Node{
        long long value, total, lazy;
        unsigned int priority;
        int size, l, r;
        bool flip;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    ImplicitTreap(const vector<long long>& values = {}){
        /// Cartesian tree on random priorities with a stack, O(n)
        vector<int> stack;
        for (long long v : values){
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

        if (!stack.empty()) root = stack[0];
        while (!stack.empty()){
            pull(stack.back());
            stack.pop_back();
        }
    }

    int new_node(long long v){
        nodes.push_back({v, v, 0, (unsigned int)rng(), 1, -1, -1, false});
        return nodes.size() - 1;
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int size() const{
        return size(root);
    }

    void apply(int t, long long delta, bool flip){
        if (t == -1) return;
        nodes[t].value += delta, nodes[t].total += delta * nodes[t].size, nodes[t].lazy += delta;
        if (flip) swap(nodes[t].l, nodes[t].r), nodes[t].flip ^= 1;
    }

    void push(int t){
        if (nodes[t].lazy || nodes[t].flip){
            apply(nodes[t].l, nodes[t].lazy, nodes[t].flip);
            apply(nodes[t].r, nodes[t].lazy, nodes[t].flip);
            nodes[t].lazy = 0, nodes[t].flip = false;
        }
    }

    void pull(int t){
        int l = nodes[t].l, r = nodes[t].r;
        nodes[t].size = 1 + size(l) + size(r);
        nodes[t].total = nodes[t].value + (l == -1 ? 0 : nodes[l].total) + (r == -1 ? 0 : nodes[r].total);
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

    void insert(int pos, long long value){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, new_node(value)), b);
    }

    void erase(int pos){
        assert(0 <= pos && pos < size());
        int a, b, c;
        split(root, pos, a, b);
        split(b, 1, b, c);
        root = merge(a, c);
    }

    /// Applies f to the treap holding exactly [l, r], then joins everything back
    template <typename F>
    void on_range(int l, int r, F f){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        f(b);
        root = merge(merge(a, b), c);
    }

    void add(int l, int r, long long delta){
        on_range(l, r, [&](int t){ apply(t, delta, false); });
    }

    void reverse(int l, int r){
        on_range(l, r, [&](int t){ apply(t, 0, true); });
    }

    long long sum(int l, int r){
        long long res = 0;
        on_range(l, r, [&](int t){ res = nodes[t].total; });
        return res;
    }

    long long get(int pos){
        return sum(pos, pos);
    }

    vector<long long> to_vector(){
        vector<long long> res;
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
};

int main(){
    Treap<int> set;
    for (int x : {5, 1, 9, 5, -1, 7}) set.insert(x);
    assert(set.size() == 6);
    assert(set.kth(0) == -1 && set.kth(1) == 1 && set.kth(2) == 5 && set.kth(3) == 5 && set.kth(5) == 9);
    assert(set.count_less(5) == 2 && set.count_less(6) == 4 && set.count_less(-5) == 0 && set.count_less(100) == 6);
    assert(set.count(5) == 2 && set.count(4) == 0);
    assert(set.erase(5) && set.count(5) == 1 && set.size() == 5);
    assert(!set.erase(4) && set.size() == 5);
    assert(set.erase(-1) && set.kth(0) == 1);

    ImplicitTreap seq({1, 2, 3, 4, 5});
    assert(seq.sum(0, 4) == 15 && seq.sum(1, 3) == 9);
    seq.add(1, 3, 10);
    assert((seq.to_vector() == vector<long long>{1, 12, 13, 14, 5}));
    seq.reverse(0, 3);
    assert((seq.to_vector() == vector<long long>{14, 13, 12, 1, 5}));
    seq.insert(2, 100);
    assert((seq.to_vector() == vector<long long>{14, 13, 100, 12, 1, 5}));
    seq.insert(6, -7);
    seq.erase(0);
    assert((seq.to_vector() == vector<long long>{13, 100, 12, 1, 5, -7}));
    assert(seq.get(1) == 100 && seq.sum(2, 5) == 11 && seq.size() == 6);

    ImplicitTreap empty;
    assert(empty.size() == 0);
    empty.insert(0, 42);
    assert(empty.get(0) == 42);

    return 0;
}
