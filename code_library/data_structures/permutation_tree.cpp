/***
 *
 * Permutation Tree
 * Splits a permutation into its common intervals, the ranges whose values form a contiguous set
 *
 * Complexity: O(n log n) to build, O(n) for count_intervals
 *
 * PermutationTree pt(p): the tree of a permutation p of 0..n - 1, n may be 0
 * pt.count_intervals(): number of common intervals, the n singletons included
 * pt.root: the root node, -1 if n = 0
 * pt.parent[v]: parent of v, -1 for the root
 * pt.children[v]: children of v in position order
 * pt.span[v]: the positions [l, r] that v covers
 * pt.range[v]: the values [lo, hi] of p[l..r], always hi - lo == r - l
 * pt.type[v]: INCREASING, DECREASING or CUT
 *
 * Notes:
 *   - Nodes 0..n - 1 are the leaves, node i is position i, internal nodes get ids n, n + 1, ..., at most 2n - 1 in all
 *   - INCREASING and DECREASING are join nodes: children ranges ascend or descend one step each
 *   - CUT is a leaf, or a node with 4+ children where only the whole node is a common interval
 *   - No join node has a child that is a join node of the same direction, which makes the tree unique
 *   - Every common interval is a leaf, a whole cut node, or a run of 2+ consecutive children of a join node
 *   - For 1-indexed or arbitrary distinct values, compress them to ranks 0..n - 1 first
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PermutationTree{
    enum Type{CUT, INCREASING, DECREASING};

    int n, root;
    vector<int> parent;
    vector<Type> type;
    vector<pair<int, int>> span, range;
    vector<vector<int>> children;

    PermutationTree(const vector<int>& p) : n(p.size()), root(-1){
        parent.reserve(2 * n), type.reserve(2 * n), span.reserve(2 * n), range.reserve(2 * n), children.reserve(2 * n);
        for (int i = 0; i < n; i++) new_node(CUT, {i, i}, {p[i], p[i]});

        /// tree_min at j holds (max - min of p[j..i]) - (i - j), which is 0 exactly when [j, i] is a common interval
        vector<int> tree_min(4 * n), tree_add(4 * n);
        vector<int> max_stack, min_stack, roots_left_to_right;
        for (int i = 0; i < n; i++){
            while (!max_stack.empty() && p[max_stack.back()] < p[i]){
                int r = max_stack.back();
                max_stack.pop_back();
                add(tree_min, tree_add, 1, 0, n - 1, max_stack.empty() ? 0 : max_stack.back() + 1, r, p[i] - p[r]);
            }
            max_stack.push_back(i);

            while (!min_stack.empty() && p[min_stack.back()] > p[i]){
                int r = min_stack.back();
                min_stack.pop_back();
                add(tree_min, tree_add, 1, 0, n - 1, min_stack.empty() ? 0 : min_stack.back() + 1, r, p[r] - p[i]);
            }
            min_stack.push_back(i);

            int cur = i;
            while (true){
                if (!roots_left_to_right.empty() && (adjacent(roots_left_to_right.back(), cur) || adjacent(cur, roots_left_to_right.back()))){
                    int top = roots_left_to_right.back();
                    Type dir = adjacent(top, cur) ? INCREASING : DECREASING;
                    roots_left_to_right.pop_back();
                    if (type[top] == dir){
                        attach(top, cur);
                        cur = top;
                        continue;
                    }

                    int v = new_node(dir, span[top], range[top]);
                    attach(v, top);
                    attach(v, cur);
                    cur = v;
                    continue;
                }

                int start = span[cur].first;
                if (start == 0 || query_min(tree_min, tree_add, 1, 0, n - 1, 0, start - 1) != 0) break;

                int v = new_node(CUT, span[cur], range[cur]);
                attach(v, cur);
                do {
                    attach(v, roots_left_to_right.back());
                    roots_left_to_right.pop_back();
                } while (range[v].second - range[v].first != span[v].second - span[v].first);
                reverse(children[v].begin(), children[v].end());
                cur = v;
            }
            roots_left_to_right.push_back(cur);

            /// Moving to i + 1 lengthens every [j, i] by one, so the -(i - j) term drops by 1 for every j <= i
            add(tree_min, tree_add, 1, 0, n - 1, 0, i, -1);
        }

        if (n) root = roots_left_to_right.back();
    }

    long long count_intervals() const {
        long long total = 0;
        for (int v = 0; v < (int)children.size(); v++){
            long long k = children[v].size();
            total += type[v] == CUT ? 1 : k * (k - 1) / 2;
        }
        return total;
    }

    int new_node(Type t, pair<int, int> s, pair<int, int> r){
        parent.push_back(-1);
        type.push_back(t);
        span.push_back(s);
        range.push_back(r);
        children.emplace_back();
        return parent.size() - 1;
    }

    static void add(vector<int>& tree_min, vector<int>& tree_add, int idx, int a, int b, int l, int r, int val){
        if (r < a || b < l) return;
        if (l <= a && b <= r){
            tree_min[idx] += val, tree_add[idx] += val;
            return;
        }

        int c = (a + b) >> 1;
        add(tree_min, tree_add, idx << 1, a, c, l, r, val);
        add(tree_min, tree_add, idx << 1 | 1, c + 1, b, l, r, val);
        tree_min[idx] = min(tree_min[idx << 1], tree_min[idx << 1 | 1]) + tree_add[idx];
    }

    /// Order matters: adjacent(left, right) means the merge of left then right is INCREASING
    bool adjacent(int u, int v) const {
        return range[u].second + 1 == range[v].first;
    }

    void attach(int v, int child){
        parent[child] = v;
        children[v].push_back(child);
        span[v] = {min(span[v].first, span[child].first), max(span[v].second, span[child].second)};
        range[v] = {min(range[v].first, range[child].first), max(range[v].second, range[child].second)};
    }

    static int query_min(const vector<int>& tree_min, const vector<int>& tree_add, int idx, int a, int b, int l, int r){
        if (r < a || b < l) return INT_MAX;
        if (l <= a && b <= r) return tree_min[idx];

        int c = (a + b) >> 1;
        return min(query_min(tree_min, tree_add, idx << 1, a, c, l, r),
                   query_min(tree_min, tree_add, idx << 1 | 1, c + 1, b, l, r)) + tree_add[idx];
    }
};

int main(){
    /// 0 [3 [1 2]] 4: an increasing root over 0, the block 3 1 2 and 4
    PermutationTree pudding({0, 3, 1, 2, 4});
    assert(pudding.count_intervals() == 10);  /// 5 singletons, 3 runs of root children, 3 1 2 and 1 2
    assert(pudding.type[pudding.root] == PermutationTree::INCREASING);
    assert(pudding.children[pudding.root].size() == 3);
    int middle = pudding.children[pudding.root][1];
    assert(pudding.type[middle] == PermutationTree::DECREASING);  /// 3 above [1 2]
    assert((pudding.span[middle] == pair<int, int>{1, 3}));

    /// No common interval besides the singletons and the whole
    PermutationTree simple({1, 3, 0, 2});
    assert(simple.type[simple.root] == PermutationTree::CUT);
    assert(simple.count_intervals() == 5);

    return 0;
}
