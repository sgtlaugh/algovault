/***
 *
 * Disjoint Set Union
 * Union-find: merges sets of nodes and tells whether two nodes share a set
 *
 * Complexity: amortized O(α(n)) per operation, union by size with path halving
 *
 * DSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they are already one set
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 * dsu.find(a): the root of a's set, the same node for every member
 *
 * Variants: rollback_dsu.cpp (undoable unions), potential_dsu.cpp (relations x[b] = x[a] + d), partially_persistent_dsu.cpp (queries at past times)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DSU{
    vector<int> sz, parent;

    DSU(int n){
        parent.resize(n + 1);
        sz.resize(n + 1, 1);
        for (int i = 0; i <= n; i++) parent[i] = i;
    }

    int find(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return parent[i];
    }

    bool connect(int a, int b){
        int c = find(a), d = find(b);
        if (c == d) return false;
        if (sz[c] > sz[d]) swap(c, d);
        parent[c] = d;
        sz[d] += sz[c], sz[c] = 0;
        return true;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        int p = find(i);
        return sz[p];
    }
};

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
    int n, m, t = 0;

    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        vector<int> w(n + 1), u(m), v(m);
        for (int i = 1; i <= n; i++) scanf("%d", &w[i]);
        for (int i = 0; i < m; i++) scanf("%d %d", &u[i], &v[i]);

        vector<array<int, 3>> ops;
        vector<char> deleted(m, 0);
        int queries = 0;
        char op;
        while (scanf(" %c", &op) == 1 && op != 'E'){
            int x, y = 0;
            if (op == 'D'){
                scanf("%d", &x);
                deleted[--x] = 1;
                ops.push_back({0, x, 0});
            }
            else if (op == 'Q'){
                scanf("%d %d", &x, &y);
                queries++;
                ops.push_back({1, x, y});
            }
            else{
                scanf("%d %d", &x, &y);
                ops.push_back({2, x, w[x]});
                w[x] = y;
            }
        }

        DSU dsu(n);
        Treap<pair<int, int>> weights;
        vector<vector<int>> members(n + 1);
        for (int i = 1; i <= n; i++) weights.insert({i, w[i]}), members[i] = {i};

        auto join = [&](int a, int b){
            int c = dsu.find(a), d = dsu.find(b);
            if (!dsu.connect(c, d)) return;
            int root = dsu.find(c), child = c ^ d ^ root;
            for (int x : members[child]){
                weights.erase({child, w[x]});
                weights.insert({root, w[x]});
                members[root].push_back(x);
            }
            members[child].clear();
        };

        for (int i = 0; i < m; i++){
            if (!deleted[i]) join(u[i], v[i]);
        }

        long long sum = 0;
        for (int i = ops.size() - 1; i >= 0; i--){
            int type = ops[i][0], x = ops[i][1], y = ops[i][2];
            if (type == 0){
                join(u[x], v[x]);
                continue;
            }

            int r = dsu.find(x);
            if (type == 1){
                int end = weights.count_less({r + 1, INT_MIN});
                if (1 <= y && y <= dsu.size(r)) sum += weights.kth(end - y).second;
            }
            else{
                weights.erase({r, w[x]});
                w[x] = y;
                weights.insert({r, w[x]});
            }
        }

        printf("Case %d: %.6f\n", ++t, (double)sum / queries);
    }
    return 0;
}
