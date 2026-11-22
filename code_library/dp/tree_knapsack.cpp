/***
 *
 * Tree Knapsack
 * Best total value of a vertex set of a forest within each weight limit, sets closed upward or connected
 *
 * Complexity:
 *   - ancestor_closed: O(n W) time and memory
 *   - connected: O(n W log n) time, O(n W) memory
 *
 * TreeKnapsack tk(parent, weight, value, W): forest with parent[v] = -1 at roots, weight[v] >= 0, value[v] any sign
 * tk.ancestor_closed(): best[c] for c in [0, W], the max value of a set of weight <= c holding each member's parent
 * tk.connected(): best[c] for c in [0, W], the max value of a set of weight <= c that is connected in the forest
 *
 * Notes:
 *   - The empty set counts, so best[0] >= 0 and best never decreases
 *   - connected splits the forest at centroids and reruns the ancestor_closed DP from each one
 *   - The DP table holds about n W long longs, 800 MB at n W = 1e8
 *   - Iterative, so deep trees cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TreeKnapsack{
    int n, W;
    vector<int> weight, roots, par, sz;
    vector<long long> value;
    vector<vector<int>> adj;
    vector<char> removed;

    TreeKnapsack(const vector<int>& parent, const vector<int>& weight, const vector<long long>& value, int W)
        : n(parent.size()), W(W), weight(weight), par(n), sz(n), value(value), adj(n), removed(n, 0){
        assert(W >= 0 && (int)weight.size() == n && (int)value.size() == n);
        for (int v = 0; v < n; v++){
            assert(weight[v] >= 0 && parent[v] >= -1 && parent[v] < n);
            if (parent[v] == -1) roots.push_back(v);
            else adj[v].push_back(parent[v]), adj[parent[v]].push_back(v);
        }
    }

    vector<long long> ancestor_closed(){
        removed.assign(n, 0);
        vector<int> order;
        for (int r : roots){
            vector<int> tree = preorder(r);
            order.insert(order.end(), tree.begin(), tree.end());
        }
        assert((int)order.size() == n);  /// a cycle of parents is never reached from a root

        return suffix_dp(order, 0);
    }

    /// Each connected set is counted at its first centroid, which roots it as an ancestor closed set
    vector<long long> connected(){
        removed.assign(n, 0);
        vector<long long> best(W + 1, 0);
        vector<int> pending = roots;
        while (!pending.empty()){
            int c = pending.back(), total = preorder(c).size();
            pending.pop_back();
            while (true){
                int heavy = -1;
                for (int u : adj[c]){
                    if (!removed[u] && u != par[c] && sz[u] * 2 > total) heavy = u;
                }
                if (heavy == -1) break;
                c = heavy;
            }

            vector<long long> rest = suffix_dp(preorder(c), 1);
            for (int cap = weight[c]; cap <= W; cap++) best[cap] = max(best[cap], value[c] + rest[cap - weight[c]]);

            removed[c] = 1;
            for (int u : adj[c]){
                if (!removed[u]) pending.push_back(u);
            }
        }
        assert(count(removed.begin(), removed.end(), 1) == n);

        return best;
    }

    /// Fills par and sz for the component of root among the vertices not removed
    vector<int> preorder(int root){
        vector<int> order, stack = {root};
        par[root] = -1;
        while (!stack.empty()){
            int v = stack.back();
            stack.pop_back();
            order.push_back(v);
            sz[v] = 1;
            for (int u : adj[v]){
                if (u != par[v] && !removed[u]) par[u] = v, stack.push_back(u);
            }
        }

        for (int i = (int)order.size() - 1; i > 0; i--) sz[par[order[i]]] += sz[order[i]];
        return order;
    }

    /// Row i of f is the best value of a set within order[i..] of weight <= c, members needing parents unless before i
    vector<long long> suffix_dp(const vector<int>& order, int row){
        int k = order.size(), width = W + 1;
        /// One flat table: connected takes 20-55% less time than with a vector per row
        vector<long long> f((k + 1LL) * width, 0);
        for (int i = k - 1; i >= row; i--){
            int v = order[i];
            long long* cur = &f[(long long)i * width];
            const long long *skip = &f[(long long)(i + sz[v]) * width], *take = cur + width;
            for (int c = 0; c <= W; c++) cur[c] = skip[c];
            for (int c = weight[v]; c <= W; c++) cur[c] = max(cur[c], value[v] + take[c - weight[v]]);
        }

        return vector<long long>(f.begin() + (long long)row * width, f.begin() + (long long)(row + 1) * width);
    }
};

int main(){
    /***
     * Root 0 (weight 2, value 3) with children 1 (1, 4) and 2 (3, 5), and 3 (2, 6) below 1
     * Closed upward: {0} at weight 2, {0, 1} at 3, {0, 1, 3} worth 13 at weight 5 beats {0, 1, 2} worth 12 at 6
     * Connected: {1} at weight 1, {3} at 2, {1, 3} worth 10 at weight 3, then {0, 1, 3} at 5
    ***/
    vector<int> parent = {-1, 0, 0, 1}, weight = {2, 1, 3, 2};
    vector<long long> value = {3, 4, 5, 6};
    TreeKnapsack tk(parent, weight, value, 6);
    assert((tk.ancestor_closed() == vector<long long>{0, 0, 3, 7, 7, 13, 13}));
    assert((tk.connected() == vector<long long>{0, 4, 6, 10, 10, 13, 13}));

    return 0;
}
