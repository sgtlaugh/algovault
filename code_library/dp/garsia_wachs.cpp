/***
 *
 * Garsia-Wachs
 * Minimum total cost of merging adjacent piles into one, a merge costing the merged size
 *
 * Complexity: expected O(n log n)
 *
 * garsia_wachs(weights): minimum total merge cost, 0 for at most one pile
 *
 * Notes:
 *   - The same value is the optimal alphabetic tree: min sum weights[i] * depth(i) over trees with the leaves in order
 *   - Weights are non-negative integers with total * max(1, ceil(log2 n)) < LLONG_MAX, which bounds the answer
 *   - Merge costs that are not the merged size: knuth_optimization.cpp
 *   - https://en.wikipedia.org/wiki/Garsia%E2%80%93Wachs_algorithm
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GarsiaWachs{
    struct Node{
        long long val, max_val;
        unsigned int priority;
        int l, r;
    };

    /// A wall above every merged pile stops the leftward walk, and end above every weight merges the last pairs
    static constexpr long long wall = LLONG_MAX, end = LLONG_MAX - 1;

    vector<Node> nodes;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};
    long long cost = 0;

    template <typename T>
    GarsiaWachs(const vector<T>& weights){
        nodes.reserve(weights.size() + 3);
        vector<int> rest = {node(end)};
        for (int i = (int)weights.size() - 1; i >= 0; i--) rest.push_back(node(weights[i]));

        /// The sequence is root then a, b; it never holds s[i] <= s[i + 2], which would merge s[i] and s[i + 1]
        int root = -1, a = node(wall), b = node(wall);
        while (!rest.empty()){
            int v = rest.back();
            rest.pop_back();

            /// Children mark a run of smaller piles put back whole; a run of one is compared as a pile and never merges
            if (nodes[v].l != -1 || nodes[v].r != -1){
                int last = pop_last(v), before_last = pop_last(v);
                root = merge(merge(merge(root, a), b), v);
                a = before_last, b = last;
            }
            else if (nodes[a].val <= nodes[v].val){
                nodes[a].val = nodes[a].max_val = nodes[a].val + nodes[b].val;
                cost += nodes[a].val;

                /// The merged pile moves left past the smaller piles, which need no new comparison except with v
                int smaller = -1;
                split_after_last_at_least(root, nodes[a].val, root, smaller);
                rest.push_back(v);
                if (smaller != -1) rest.push_back(smaller);
                rest.push_back(a);
                b = pop_last(root), a = pop_last(root);
            }
            else{
                root = merge(root, a);
                a = b, b = v;
            }
        }
    }

    int node(long long x){
        nodes.push_back({x, x, (unsigned int)rng(), -1, -1});
        return nodes.size() - 1;
    }

    int merge(int x, int y){
        if (x == -1 || y == -1) return x == -1 ? y : x;
        if (nodes[x].priority > nodes[y].priority){
            nodes[x].r = merge(nodes[x].r, y);
            pull(x);
            return x;
        }
        nodes[y].l = merge(x, nodes[y].l);
        pull(y);
        return y;
    }

    /// The returned node has no children, or the main loop would read it as a run of piles
    int pop_last(int& t){
        if (nodes[t].r == -1){
            int last = t;
            t = nodes[t].l;
            nodes[last].l = -1, nodes[last].max_val = nodes[last].val;
            return last;
        }
        int last = pop_last(nodes[t].r);
        pull(t);
        return last;
    }

    void pull(int t){
        nodes[t].max_val = nodes[t].val;
        if (nodes[t].l != -1) nodes[t].max_val = max(nodes[t].max_val, nodes[nodes[t].l].max_val);
        if (nodes[t].r != -1) nodes[t].max_val = max(nodes[t].max_val, nodes[nodes[t].r].max_val);
    }

    void split_after_last_at_least(int t, long long at_least, int& x, int& y){
        if (t == -1){
            x = y = -1;
            return;
        }

        if (nodes[t].r != -1 && nodes[nodes[t].r].max_val >= at_least){
            split_after_last_at_least(nodes[t].r, at_least, nodes[t].r, y);
            x = t;
        }
        else if (nodes[t].val >= at_least){
            y = nodes[t].r, nodes[t].r = -1;
            x = t;
        }
        else{
            split_after_last_at_least(nodes[t].l, at_least, x, nodes[t].l);
            y = t;
        }
        pull(t);
    }
};

template <typename T>
long long garsia_wachs(const vector<T>& weights){
    return GarsiaWachs(weights).cost;
}

int main(){
    /// 2 + 3 = 5, then 5 + 4 = 9, then 5 + 9 = 14
    assert(garsia_wachs(vector<int>{5, 2, 3, 4}) == 28);
    assert(garsia_wachs(vector<int>{1, 2, 3}) == 9);
    assert(garsia_wachs(vector<int>{4, 1, 1, 4}) == 18);
    assert(garsia_wachs(vector<long long>(8, 1)) == 24);

    assert(garsia_wachs(vector<int>{7}) == 0);
    assert(garsia_wachs(vector<int>{}) == 0);
    assert(garsia_wachs(vector<int>{0, 0, 0}) == 0);

    return 0;
}
