/***
 *
 * Sprague-Grundy
 * Grundy values of impartial games, a position is lost for the player to move iff its value is 0
 *
 * Complexity:
 *   - Mex: O(1) per add and clear, O(answer) per value
 *   - grundy_dag: O(n + m) for n vertices and m edges
 *   - grundy_table: O(n + total moves)
 *
 * Mex mex: an empty set of non-negative integers
 * mex.add(x): inserts x
 * mex.clear(): empties the set
 * mex.value(): the smallest non-negative integer not in the set
 *
 * grundy_dag(adj): g[u] for each vertex u of a DAG, the mex of g[v] over the edges u -> v in adj[u]
 * grundy_table(n, moves): g[0..n], g[i] the Grundy value of a heap of size i, the mex of what moves(i, g, mex) adds
 *
 * Notes:
 *   - moves(i, g, mex) calls mex.add(x) once per move from a heap of size i, x the value it leads to, g[0..i - 1] known
 *   - Independent parts combine by xor: a move splitting a heap into sizes a and b adds g[a] ^ g[b]
 *   - A position of several heaps or tokens is worth the xor of their values
 *   - Values pass 63: Treblecross reaches 262 below n = 1e4, so a 64-bit seen mask overflows silently
 *   - Mex uses memory linear in the largest x added, x >= 0
 *   - grundy_dag asserts that adj has no cycle, self loops and parallel edges included
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// seen[x] == stamp marks x, so clear is O(1): 1.5x faster than filling a vector and scanning it, n = 1e4 Treblecross
struct Mex{
    vector<int> seen;
    int stamp = 1;

    void add(int x){
        if (x >= (int)seen.size()) seen.resize(2 * x + 2);
        seen[x] = stamp;
    }

    void clear(){
        stamp++;
    }

    int value() const{
        int res = 0;
        while (res < (int)seen.size() && seen[res] == stamp) res++;
        return res;
    }
};

vector<int> grundy_dag(const vector<vector<int>>& adj){
    int n = adj.size();
    vector<vector<int>> parents(n);
    vector<int> left(n), order;
    for (int u = 0; u < n; u++){
        left[u] = adj[u].size();
        if (!left[u]) order.push_back(u);
        for (int v : adj[u]) parents[v].push_back(u);
    }

    for (int i = 0; i < (int)order.size(); i++){
        for (int p : parents[order[i]]){
            if (--left[p] == 0) order.push_back(p);
        }
    }
    assert((int)order.size() == n);

    vector<int> g(n);
    Mex mex;
    for (int u : order){
        mex.clear();
        for (int v : adj[u]) mex.add(g[v]);
        g[u] = mex.value();
    }

    return g;
}

template <typename Moves>
vector<int> grundy_table(int n, Moves moves){
    vector<int> g(n + 1);
    Mex mex;
    for (int i = 0; i <= n; i++){
        mex.clear();
        moves(i, g, mex);
        g[i] = mex.value();
    }

    return g;
}

int main(){
    /// Subtraction game, take 1 or 2 from the heap: the values cycle 0, 1, 2
    vector<int> take = grundy_table(7, [](int i, const vector<int>& g, Mex& mex){
        if (i >= 1) mex.add(g[i - 1]);
        if (i >= 2) mex.add(g[i - 2]);
    });
    assert((take == vector<int>{0, 1, 2, 0, 1, 2, 0, 1}));

    /// Kayles: knock down one pin or two adjacent pins of a row, splitting it in two rows
    vector<int> kayles = grundy_table(5, [](int i, const vector<int>& g, Mex& mex){
        for (int j = 0; j + 1 <= i; j++) mex.add(g[j] ^ g[i - 1 - j]);
        for (int j = 0; j + 2 <= i; j++) mex.add(g[j] ^ g[i - 2 - j]);
    });
    assert((kayles == vector<int>{0, 1, 2, 3, 1, 4}));
    assert((kayles[4] ^ kayles[5] ^ kayles[1]) == 4);  /// rows of 4, 5 and 1 pins: the first player wins

    /// Tokens on a DAG, a move slides one token along an edge
    vector<vector<int>> adj = {{1, 2}, {2}, {3}, {}, {0, 3}};
    vector<int> g = grundy_dag(adj);
    assert((g == vector<int>{2, 0, 1, 0, 1}));
    assert((g[0] ^ g[4]) == 3);

    /// A value far above 63 is stored, not wrapped into a 64-bit mask
    Mex mex;
    mex.add(1);
    mex.add(1000);
    mex.add(0);
    assert(mex.value() == 2);

    mex.add(2);
    assert(mex.value() == 3);

    mex.clear();
    mex.add(64);
    assert(mex.value() == 0);

    return 0;
}
