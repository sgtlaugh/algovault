/***
 *
 * Hackenbush
 * Green Hackenbush Grundy value of a rooted graph and red-blue Hackenbush value of a rooted tree
 *
 * Complexity:
 *   - GreenHackenbush grundy: O(n + m) per call
 *   - RedBlueHackenbush value: O(n log^2 n + m) per call
 *
 * GreenHackenbush g(n): impartial game on nodes 0 to n - 1, either player deletes any edge
 * g.add_edge(u, v): parallel edges and self loops (u == v) allowed
 * g.grundy(root): the Grundy value, the first player wins iff it is nonzero
 *
 * RedBlueHackenbush t(n): partizan game on nodes 0 to n - 1, Left deletes only blue edges, Right only red
 * t.add_edge(u, v, blue): adds a blue edge if blue, else a red one
 * t.value(root): the exact value as a Dyadic, the root's component must be a tree
 *
 * Dyadic: whole + sum of 2^-b over b in bits, whole = floor, bits strictly increasing, all >= 1
 * d.sign(): > 0 if Left wins whoever starts, < 0 if Right wins whoever starts, 0 if the second player wins
 *
 * Notes:
 *   - A move deletes one edge and every edge it cuts off from the root, a player who cannot move loses
 *   - Several ground vertices: map them all to one root vertex before add_edge
 *   - Edges outside the root's component take no part in the game
 *   - In RedBlueHackenbush parallel edges and self loops are cycles
 *   - Denominators reach 2^(n - 2) on a long stalk, so the value is kept as bits, not as a double
 *   - Fusion principle: the vertices of a 2-edge-connected part fuse into one, its edges become loops (1 each)
 *   - Colon principle: a subtree of value x on top of a bridge is worth x + 1, siblings combine by xor
 *   - A blue edge over a subtree of value x is worth x + 1 for x >= 0, else (1 + frac(x)) / 2^(-floor(x))
 *   - A red edge over x is worth the negation of the blue value of -x, and sibling values add
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GreenHackenbush{
    struct Frame{
        int u, up, next, end, dis, low, value;
    };

    int n;
    vector<pair<int, int>> edges;

    GreenHackenbush(int n) : n(n) {}

    void add_edge(int u, int v){
        edges.push_back({u, v});
    }

    int grundy(int root){
        vector<int> degree(n, 0), neighbors(n, 0), value(n, 0); /// neighbors xors the ids, the only one at degree 1
        for (auto [u, v] : edges){
            if (u == v){
                value[u] ^= 1;
                continue;
            }
            degree[u]++, neighbors[u] ^= v;
            degree[v]++, neighbors[v] ^= u;
        }

        /// Pendant trees fold into the vertex they hang from by the colon principle, with no adjacency lists built
        vector<int> leaves;
        for (int u = 0; u < n; u++){
            if (u != root && degree[u] == 1) leaves.push_back(u);
        }
        for (int i = 0; i < (int)leaves.size(); i++){
            int u = leaves[i], p = neighbors[u];
            if (degree[u] != 1) continue; /// its last neighbor was folded into it, in a component without the root

            degree[u] = 0;
            degree[p]--, neighbors[p] ^= u, value[p] ^= value[u] + 1;
            if (p != root && degree[p] == 1) leaves.push_back(p);
        }
        if (degree[root] == 0) return value[root];

        /// What is left is cycles and paths between them: Tarjan bridges on the edges with both ends unfolded
        vector<int> start(n + 1, 0), adj, dis(n, -1);
        for (auto [u, v] : edges){
            if (degree[u] && degree[v]) start[u]++, start[v]++;
        }
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        adj.resize(start[n]);
        for (auto [u, v] : edges){
            if (degree[u] && degree[v]) adj[--start[u]] = v, adj[--start[v]] = u;
        }

        int timer = 0;
        vector<Frame> stack = {{root, -1, start[root], start[root + 1], timer, timer, value[root]}};
        dis[root] = timer++;

        while (true){
            Frame& f = stack.back();
            if (f.next < f.end){
                int v = adj[f.next++];
                /// One edge to the parent is the tree edge, any parallel copy of it closes a cycle
                if (v == f.up) f.up = -1;
                else if (dis[v] == -1){
                    dis[v] = timer;
                    stack.push_back({v, f.u, start[v], start[v + 1], timer, timer, value[v]});
                    timer++;
                }
                else if (dis[v] < f.dis){
                    f.low = min(f.low, dis[v]);
                    f.value ^= 1; /// a cycle edge fuses into a loop worth 1, counted once, at its deeper endpoint
                }
                continue;
            }

            if (stack.size() == 1) return f.value;
            Frame& p = stack[stack.size() - 2];
            p.low = min(p.low, f.low);
            p.value ^= f.low > p.dis ? f.value + 1 : f.value ^ 1;
            stack.pop_back();
        }
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

        /// (frac - 2) / 2^(w + 1) = -1 + (1 - 2^-w) + frac / 2^(w + 1): w new bits, paid by the w integer part edges
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

    int n;
    vector<array<int, 3>> edges; /// (u, v, blue)

    RedBlueHackenbush(int n) : n(n) {}

    void add_edge(int u, int v, bool blue){
        edges.push_back({u, v, blue});
    }

    Dyadic value(int root){
        vector<int> degree(n, 0), incident(n, 0); /// incident xors the edge ids, the only one at degree 1
        for (int id = 0; id < (int)edges.size(); id++){
            auto [u, v, blue] = edges[id];
            degree[u]++, incident[u] ^= id;
            degree[v]++, incident[v] ^= id;
        }

        /// Leaves fold into their parent until only the root is left, so no adjacency lists are built
        vector<int> leaves;
        for (int u = 0; u < n; u++){
            if (u != root && degree[u] == 1) leaves.push_back(u);
        }

        vector<Number> num(n);
        for (int i = 0; i < (int)leaves.size(); i++){
            int u = leaves[i];
            if (degree[u] != 1) continue; /// its last neighbor was folded into it, in a component without the root

            auto [a, b, blue] = edges[incident[u]];
            int p = a ^ b ^ u;
            degree[u] = 0;
            degree[p]--, incident[p] ^= incident[u];
            if (p != root && degree[p] == 1) leaves.push_back(p);

            blue ? num[u].below_blue() : num[u].below_red();
            if (num[p].bits.size() < num[u].bits.size()) swap(num[p], num[u]);
            num[p].absorb(num[u]);
            num[u] = Number();
        }
        assert(degree[root] == 0); /// a cycle in the root's component never folds away, so the root keeps an edge

        Dyadic result;
        result.whole = num[root].whole;
        for (int b : num[root].bits) result.bits.push_back(b + num[root].shift);
        return result;
    }
};

int main(){
    /// A stalk of k edges is a nim heap of size k, siblings xor
    GreenHackenbush path(5);
    for (int i = 0; i + 1 < 5; i++) path.add_edge(i, i + 1);
    assert(path.grundy(0) == 4);
    assert(path.grundy(2) == 0);                /// grounded in the middle: 2 ^ 2

    /// Fusion: a cycle of k edges becomes k loops, worth k & 1
    GreenHackenbush triangle(3);
    for (int i = 0; i < 3; i++) triangle.add_edge(i, (i + 1) % 3);
    assert(triangle.grundy(0) == 1);

    /// Edge 0-1 holding the triangle 1 2 3, the colon principle gives 1 + 1
    GreenHackenbush lollipop(4);
    lollipop.add_edge(0, 1);
    for (int i = 1; i <= 3; i++) lollipop.add_edge(i, i % 3 + 1);
    assert(lollipop.grundy(0) == 2);

    auto stalk = [](const string& colors){
        RedBlueHackenbush t(colors.size() + 1);
        for (int i = 0; i < (int)colors.size(); i++) t.add_edge(i, i + 1, colors[i] == 'B');
        return t.value(0);
    };
    assert((stalk("BB") == Dyadic{2, {}}));
    assert((stalk("BR") == Dyadic{0, {1}}));    /// 1/2
    assert((stalk("BRB") == Dyadic{0, {1, 2}}));  /// 1/2 + 1/4
    assert(stalk("RB").sign() == -1);           /// -1/2, Right wins whoever starts

    /// A blue-red stalk beside a red-blue stalk: 1/2 - 1/2, the second player wins
    RedBlueHackenbush balanced(5);
    balanced.add_edge(0, 1, true);
    balanced.add_edge(1, 2, false);
    balanced.add_edge(0, 3, false);
    balanced.add_edge(3, 4, true);
    assert(balanced.value(0).sign() == 0);
    return 0;
}
