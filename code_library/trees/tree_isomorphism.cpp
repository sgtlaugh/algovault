/***
 *
 * Tree Isomorphism (AHU canonical ids)
 * Rooted and unrooted isomorphism of trees, 0-based nodes
 *
 * Complexity: O(n log N) per tree, N = total nodes canonized by the instance, O(N) memory for the dictionary
 *
 * Every rooted subtree gets an integer id from a dictionary keyed by the sorted ids of its children,
 * so two rooted trees are isomorphic iff their roots get the same id from the same TreeCanonizer
 * An unrooted tree's id is the smaller rooted id over its 1 or 2 centers, since isomorphisms map centers to centers
 * The dictionary is a map, so equal ids are exact isomorphism, never a hash collision
 *
 * Usage:
 *   IsoTree a(n), b(n); a.add_edge(u, v); ...
 *   TreeCanonizer canon;
 *   canon.rooted_id(a, ra) == canon.rooted_id(b, rb)    rooted isomorphism
 *   canon.unrooted_id(a) == canon.unrooted_id(b)         unrooted isomorphism
 *   canon.subtree_ids(a, root)[v]                        id of the subtree of v, distinct values = distinct subtrees
 *   a.centers()                                          the 1 or 2 centers, ascending
 *
 * Ids are only comparable within one TreeCanonizer, and unrooted ids only with unrooted ids
 * The dictionary keeps every distinct subtree it has seen, use a fresh instance to free it
 * The graph must be a tree with n >= 1, every call asserts n - 1 edges and that every node is reached exactly once
 * Traversals are iterative so deep trees are safe
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
};

struct TreeCanonizer{
    map<vector<int>, int> ids;

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

    int rooted_id(const IsoTree& t, int root){
        return subtree_ids(t, root)[root];
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
    assert(canon.rooted_id(t1, 1) != canon.rooted_id(t2, 0));
    assert(canon.rooted_id(t1, 1) == canon.rooted_id(t2, 4));
    assert((t1.centers() == vector<int>{0, 1}));
    assert((t2.centers() == vector<int>{0, 4}));

    vector<int> sub = canon.subtree_ids(t1, 0);
    assert(sub[2] == sub[3] && sub[3] == sub[4]);
    assert(sub[1] != sub[2] && sub[0] != sub[1] && sub[0] != sub[2]);
    assert(set<int>(sub.begin(), sub.end()).size() == 3);

    IsoTree path = make_tree(4, {{0, 1}, {1, 2}, {2, 3}});
    IsoTree path_relabeled = make_tree(4, {{2, 0}, {0, 3}, {3, 1}});
    IsoTree star = make_tree(4, {{0, 1}, {0, 2}, {0, 3}});
    assert(canon.unrooted_id(path) == canon.unrooted_id(path_relabeled));
    assert(canon.unrooted_id(path) != canon.unrooted_id(star));
    assert((path.centers() == vector<int>{1, 2}));
    assert((path_relabeled.centers() == vector<int>{0, 3}));
    assert((star.centers() == vector<int>{0}));

    vector<int> rooted4 = {canon.rooted_id(path, 0), canon.rooted_id(path, 1), canon.rooted_id(star, 0), canon.rooted_id(star, 1)};
    assert(set<int>(rooted4.begin(), rooted4.end()).size() == 4);
    assert(canon.rooted_id(path, 3) == rooted4[0] && canon.rooted_id(path, 2) == rooted4[1]);
    assert(canon.rooted_id(star, 3) == rooted4[3]);

    /// the 6 unlabeled trees on 6 nodes (OEIS A000055)
    vector<IsoTree> six = {
        make_tree(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}}),
        make_tree(6, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}),
        make_tree(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {1, 5}}),
        make_tree(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 5}}),
        make_tree(6, {{0, 1}, {1, 2}, {2, 3}, {1, 4}, {1, 5}}),
        make_tree(6, {{0, 1}, {1, 2}, {2, 3}, {1, 4}, {2, 5}}),
    };
    set<int> six_ids;
    for (auto& t : six) six_ids.insert(canon.unrooted_id(t));
    assert(six_ids.size() == 6);

    IsoTree h_relabeled = make_tree(6, {{5, 3}, {3, 0}, {0, 2}, {3, 1}, {0, 4}});
    assert(canon.unrooted_id(h_relabeled) == canon.unrooted_id(six[5]));
    assert(canon.unrooted_id(h_relabeled) != canon.unrooted_id(six[2]));
    assert((six[0].centers() == vector<int>{2, 3}));
    assert((six[3].centers() == vector<int>{2}));

    IsoTree single(1), other_single(1);
    assert(canon.unrooted_id(single) == canon.unrooted_id(other_single));
    assert(canon.rooted_id(single, 0) == sub[2]);
    assert((single.centers() == vector<int>{0}));

    IsoTree edge = make_tree(2, {{1, 0}});
    assert((edge.centers() == vector<int>{0, 1}));
    assert(canon.rooted_id(edge, 0) == canon.rooted_id(edge, 1));
    assert(canon.unrooted_id(edge) != canon.unrooted_id(single));

    return 0;
}
