/***
 *
 * Hackenbush
 * Green Hackenbush Grundy value of a rooted graph and red-blue Hackenbush value of a rooted tree
 *
 * Both games: a move deletes one edge, then every edge no longer connected to the ground (root) disappears
 * The player who cannot move loses. Several ground vertices: map them all to one root vertex before add_edge
 * Edges outside the root's component take no part in the game. Nodes are numbered from 0 to n - 1
 * Traversals are iterative, so deep graphs are safe
 *
 * GreenHackenbush (impartial, either player deletes any edge):
 *   - Complexity: O(n + m) per grundy call
 *   - GreenHackenbush g(n); g.add_edge(u, v); g.grundy(root): first player wins iff the result is nonzero
 *   - Parallel edges and self loops (u == v) allowed
 *   - Fusion principle: the vertices of a 2-edge-connected part fuse into one, its edges become loops (1 each)
 *   - Colon principle: a subtree of value x on top of a bridge is worth x + 1, siblings combine by xor
 *
 * RedBlueHackenbush (partizan, Left deletes only blue edges, Right only red):
 *   - Complexity: O(n log^2 n) per value call, small-to-large merging of fraction bits into the heavy child
 *   - RedBlueHackenbush t(n); t.add_edge(u, v, blue); t.value(root): the exact value as a Dyadic
 *   - The root's component must be a tree (asserted), parallel edges and self loops are cycles
 *   - Dyadic: whole + sum of 2^-b over b in bits, whole = floor, bits strictly increasing, all >= 1
 *     Denominators reach 2^(n - 2) (a long stalk), so the value is kept as bits, not as a double
 *   - sign() > 0: Left wins whoever starts, < 0: Right wins whoever starts, 0: the second player wins
 *   - An edge over a subtree of value x is worth, if blue: x + 1 for x >= 0, else (1 + frac(x)) / 2^(-floor(x))
 *     If red: the negation of the blue value of -x. Siblings add
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GreenHackenbush{
    int n, m = 0;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge index)
    vector<int> loops;

    GreenHackenbush(int n) : n(n), adj(n), loops(n) {}

    void add_edge(int u, int v){
        if (u == v){
            loops[u]++;
            return;
        }

        adj[u].push_back({v, m});
        adj[v].push_back({u, m});
        m++;
    }

    int grundy(int root){
        int timer = 0;
        vector<int> dis(n, -1), low(n), parent_edge(n, -1), next(n, 0), value(n, 0), stack = {root};
        dis[root] = low[root] = timer++;
        value[root] = loops[root] & 1;

        while (!stack.empty()){
            int u = stack.back();
            if (next[u] < (int)adj[u].size()){
                auto [v, id] = adj[u][next[u]++];
                if (id == parent_edge[u]) continue;
                if (dis[v] == -1){
                    dis[v] = low[v] = timer++;
                    parent_edge[v] = id, value[v] = loops[v] & 1;
                    stack.push_back(v);
                }
                else{
                    low[u] = min(low[u], dis[v]);
                    /// A cycle edge fuses into a loop worth 1, counted once, at its deeper endpoint
                    if (dis[v] < dis[u]) value[u] ^= 1;
                }
                continue;
            }

            stack.pop_back();
            if (stack.empty()) break;
            int p = stack.back();
            low[p] = min(low[p], low[u]);
            value[p] ^= low[u] > dis[p] ? value[u] + 1 : value[u] ^ 1;
        }

        return value[root];
    }
};

struct Dyadic{
    long long whole = 0;
    vector<int> bits;

    int sign() const{
        if (whole != 0) return whole > 0 ? 1 : -1;
        return !bits.empty();
    }

    bool operator==(const Dyadic& other) const{
        return whole == other.whole && bits == other.bits;
    }
};

struct RedBlueHackenbush{
    struct Number{
        long long whole = 0;
        int shift = 0;
        set<int> bits; /// actual exponent = stored + shift, so halving every bit at once is O(1)

        void add_bit(int e){
            while (e > 0 && bits.erase(e - shift)) e--;
            if (e == 0) whole++;
            else bits.insert(e - shift);
        }

        void absorb(const Number& other){
            whole += other.whole;
            for (int b : other.bits) add_bit(b + other.shift);
        }

        void below_blue(){
            if (whole >= 0){
                whole++;
                return;
            }

            int k = -whole;
            shift += k, whole = 0;
            bits.insert(k - shift);
        }

        /// (frac - 2) / 2^(w + 1) = -1 + (1 - 2^-w) + frac / 2^(w + 1): the w new bits are paid for by the
        /// w edges that built the integer part, so all calls together insert O(n) bits
        void below_red(){
            if (whole < 0 || (whole == 0 && bits.empty())){
                whole--;
                return;
            }

            int w = whole;
            shift += w + 1, whole = -1;
            for (int e = 1; e <= w; e++) bits.insert(e - shift);
        }
    };

    int n, m = 0;
    vector<vector<array<int, 3>>> adj; /// (neighbor, edge index, blue)

    RedBlueHackenbush(int n) : n(n), adj(n) {}

    void add_edge(int u, int v, bool blue){
        adj[u].push_back({v, m, blue});
        adj[v].push_back({u, m, blue});
        m++;
    }

    Dyadic value(int root){
        vector<int> order = {root}, parent(n, -1), parent_edge(n, -1), size(n, 0);
        vector<char> visited(n, 0), blue(n, 0);
        visited[root] = 1;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (auto [v, id, b] : adj[u]){
                if (id == parent_edge[u]) continue;
                assert(!visited[v]);
                visited[v] = 1, parent[v] = u, parent_edge[v] = id, blue[v] = b;
                order.push_back(v);
            }
        }

        for (int i = (int)order.size() - 1; i > 0; i--) size[parent[order[i]]] += size[order[i]] + 1;

        vector<Number> num(n);
        for (int i = (int)order.size() - 1; i >= 0; i--){
            int u = order[i], heavy = -1;
            for (auto [v, id, b] : adj[u]){
                if (id != parent_edge[u] && (heavy == -1 || size[v] > size[heavy])) heavy = v;
            }
            if (heavy == -1) continue;

            num[u] = move(num[heavy]);
            blue[heavy] ? num[u].below_blue() : num[u].below_red();
            for (auto [v, id, b] : adj[u]){
                if (id == parent_edge[u] || v == heavy) continue;
                b ? num[v].below_blue() : num[v].below_red();
                num[u].absorb(num[v]);
                num[v] = Number();
            }
        }

        Dyadic result;
        result.whole = num[root].whole;
        for (int b : num[root].bits) result.bits.push_back(b + num[root].shift);
        return result;
    }
};

int main(){
    GreenHackenbush lonely(1);
    assert(lonely.grundy(0) == 0);

    GreenHackenbush path(5);
    for (int i = 0; i + 1 < 5; i++) path.add_edge(i, i + 1);
    assert(path.grundy(0) == 4);
    assert(path.grundy(1) == (1 ^ 3));
    assert(path.grundy(2) == 0);

    GreenHackenbush triangle(3);
    triangle.add_edge(0, 1);
    triangle.add_edge(1, 2);
    triangle.add_edge(2, 0);
    assert(triangle.grundy(0) == 1);

    GreenHackenbush square(4);
    for (int i = 0; i < 4; i++) square.add_edge(i, (i + 1) % 4);
    assert(square.grundy(0) == 0);

    GreenHackenbush loops(5);
    loops.add_edge(0, 0);
    loops.add_edge(0, 1);
    loops.add_edge(1, 1);
    loops.add_edge(3, 4);
    loops.add_edge(4, 4);
    assert(loops.grundy(0) == (1 ^ 2));
    loops.add_edge(0, 0);
    assert(loops.grundy(0) == 2);

    GreenHackenbush doubled(2);
    doubled.add_edge(0, 1);
    doubled.add_edge(1, 0);
    assert(doubled.grundy(0) == 0);

    GreenHackenbush lollipop(4);
    lollipop.add_edge(0, 1);
    lollipop.add_edge(1, 2);
    lollipop.add_edge(2, 3);
    lollipop.add_edge(3, 1);
    assert(lollipop.grundy(0) == 2);
    assert(lollipop.grundy(1) == (1 ^ 1));

    RedBlueHackenbush empty(1);
    assert(empty.value(0).sign() == 0 && empty.value(0) == Dyadic());

    auto stalk = [](const string& colors){
        RedBlueHackenbush t(colors.size() + 1);
        for (int i = 0; i < (int)colors.size(); i++) t.add_edge(i, i + 1, colors[i] == 'B');
        return t.value(0);
    };

    assert((stalk("B") == Dyadic{1, {}}));
    assert((stalk("R") == Dyadic{-1, {}}));
    assert((stalk("BR") == Dyadic{0, {1}}));
    assert((stalk("BRB") == Dyadic{0, {1, 2}}));
    assert((stalk("BRR") == Dyadic{0, {2}}));
    assert((stalk("BBR") == Dyadic{1, {1}}));
    assert((stalk("RB") == Dyadic{-1, {1}}));
    assert((stalk("RBB") == Dyadic{-1, {1, 2}}));
    assert(stalk("RB").sign() == -1 && stalk("BRR").sign() == 1);

    RedBlueHackenbush balanced(5);
    balanced.add_edge(0, 1, true);
    balanced.add_edge(1, 2, false);
    balanced.add_edge(0, 3, false);
    balanced.add_edge(3, 4, true);
    assert(balanced.value(0).sign() == 0 && balanced.value(0) == Dyadic());

    RedBlueHackenbush fork(7);
    fork.add_edge(0, 1, true);
    fork.add_edge(1, 2, false);
    fork.add_edge(1, 3, false);
    assert((fork.value(0) == Dyadic{0, {2}}));
    fork.add_edge(0, 4, true);
    fork.add_edge(4, 5, false);
    fork.add_edge(5, 6, true);
    assert((fork.value(0) == Dyadic{1, {}}));

    return 0;
}
