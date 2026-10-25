/***
 *
 * Tree Isomorphism (AHU canonical ids)
 * Rooted and unrooted isomorphism of trees, 0-based nodes
 *
 * Complexity: O(n log N) per tree, N = total nodes of all trees given to the instance, O(N) memory for the dictionary
 *
 * IsoTree a(n): an empty tree on n >= 1 nodes
 * a.add_edge(u, v): adds the edge u - v
 * a.centers(): the 1 or 2 centers, ascending
 *
 * TreeCanonizer canon: an empty dictionary of subtree ids
 * canon.rooted_id(a, root): the id of a rooted at root, equal for two rooted trees iff they are isomorphic
 * canon.subtree_ids(a, root): ids where ids[v] is the id of the subtree of v, equal iff the subtrees are isomorphic
 * canon.unrooted_id(a): the id of a as an unrooted tree, equal for two trees iff they are isomorphic
 *
 * Notes:
 *   - Ids are only comparable within one TreeCanonizer, and unrooted ids only with unrooted ids
 *   - Equal ids are exact isomorphism, never a hash collision
 *   - The dictionary keeps every distinct subtree it has seen, use a fresh instance to free it
 *   - Every call asserts n - 1 edges and that every node is reached exactly once
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct IsoTree{
    int n, added_edges = 0;
    vector<vector<int>> adj;

    IsoTree(int n) : n(n), adj(n){
        assert(n >= 1);
    }

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
        added_edges++;
    }

    /// Peels leaves layer by layer, the last layer of 1 or 2 nodes is the center
    vector<int> centers() const{
        bfs(0);
        vector<int> deg(n), layer;
        for (int v = 0; v < n; v++){
            deg[v] = adj[v].size();
            if (deg[v] <= 1) layer.push_back(v);
        }

        for (int remaining = n; remaining > 2;){
            remaining -= layer.size();
            vector<int> next;
            for (int u : layer){
                for (int v : adj[u]){
                    if (--deg[v] == 1) next.push_back(v);
                }
            }
            layer = next;
        }

        sort(layer.begin(), layer.end());
        return layer;
    }

    /// parent[root] = root
    pair<vector<int>, vector<int>> bfs(int root) const{
        assert(0 <= root && root < n);
        assert(added_edges == n - 1);
        vector<int> order = {root}, parent(n, -1);
        parent[root] = root;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int v : adj[u]){
                if (v == parent[u]) continue;
                assert(parent[v] == -1);
                parent[v] = u;
                order.push_back(v);
            }
        }

        assert((int)order.size() == n);
        return {order, parent};
    }
};

struct TreeCanonizer{
    map<vector<int>, int> ids;

    int rooted_id(const IsoTree& t, int root){
        return subtree_ids(t, root)[root];
    }

    /// Reverse BFS order finishes every child before its parent
    vector<int> subtree_ids(const IsoTree& t, int root){
        auto [order, parent] = t.bfs(root);
        vector<vector<int>> children(t.n);
        vector<int> id(t.n);

        for (int i = t.n - 1; i >= 0; i--){
            int v = order[i];
            sort(children[v].begin(), children[v].end());
            int next_id = ids.size();
            id[v] = ids.try_emplace(move(children[v]), next_id).first->second;
            if (v != root) children[parent[v]].push_back(id[v]);
        }

        return id;
    }

    int unrooted_id(const IsoTree& t){
        int res = INT_MAX;
        for (int c : t.centers()) res = min(res, rooted_id(t, c));
        return res;
    }
};

int main(){
    auto make_tree = [](int n, const vector<pair<int, int>>& edges){
        IsoTree t(n);
        for (auto [u, v] : edges) t.add_edge(u, v);
        return t;
    };

    /***
     *      t1            t2
     *       0             0
     *      / \           / \
     *     1   4         1   4
     *    / \               / \
     *   2   3             3   2
    ***/
    IsoTree t1 = make_tree(5, {{0, 1}, {1, 2}, {1, 3}, {0, 4}});
    IsoTree t2 = make_tree(5, {{0, 1}, {0, 4}, {4, 3}, {4, 2}});
    TreeCanonizer canon;

    assert(canon.unrooted_id(t1) == canon.unrooted_id(t2));
    assert(canon.rooted_id(t1, 0) == canon.rooted_id(t2, 0));
    assert(canon.rooted_id(t1, 1) == canon.rooted_id(t2, 4));   /// two leaves and a path of 2 below each
    assert(canon.rooted_id(t1, 1) != canon.rooted_id(t2, 0));
    assert((t1.centers() == vector<int>{0, 1}));

    vector<int> sub = canon.subtree_ids(t1, 0);
    assert(sub[2] == sub[4]);                                   /// every leaf is the same subtree
    assert(sub[1] != sub[2]);

    IsoTree path = make_tree(4, {{0, 1}, {1, 2}, {2, 3}});
    IsoTree star = make_tree(4, {{0, 1}, {0, 2}, {0, 3}});
    assert(canon.unrooted_id(path) != canon.unrooted_id(star));
    assert((path.centers() == vector<int>{1, 2}));
    assert((star.centers() == vector<int>{0}));
    return 0;
}
