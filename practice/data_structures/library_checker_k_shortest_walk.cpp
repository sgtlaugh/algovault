/***
 *
 * Persistent Leftist Heap
 * Meldable priority queue where every operation returns a new version and old versions stay valid
 *
 * Complexity: O(log n) time and at most 2 * log2(n + 1) + 2 new nodes per push, pop and meld, n = size of the largest heap involved
 *
 * LeftistHeap<T> heap: a pool of min-heaps, each version is an int and EMPTY (-1) is the empty heap
 * LeftistHeap<T, greater<T>> heap: the same pool of max-heaps
 * heap.push(h, x): a new version holding the values of h plus x
 * heap.pop(h): a new version holding the values of h without its top, h != EMPTY
 * heap.top(h): first value of h under Compare, the smallest for less<T>, h != EMPTY
 * heap.meld(a, b): a new version holding the values of a and b, meld(h, h) holds every value twice
 *
 * No call changes an existing version, and a version only works with the object that made it
 * Nodes are never freed, so q operations hold O(q log n) nodes
 * Usual use: k shortest walks (Eppstein), where each vertex's heap extends its parent's
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED LeftistHeap
template <typename T, typename Compare = less<T>>
struct LeftistHeap{
    static constexpr int EMPTY = -1;

    struct Node{
        T val;
        int left, right, rank;
    };

    vector<Node> nodes;
    Compare cmp;

    int push(int h, const T& x){
        nodes.push_back({x, EMPTY, EMPTY, 1});
        return meld(h, (int)nodes.size() - 1);
    }

    int pop(int h){
        assert(h != EMPTY);
        return meld(nodes[h].left, nodes[h].right);
    }

    T top(int h) const{
        assert(h != EMPTY);
        return nodes[h].val;
    }

    /// Copies only the nodes on the merged right spines, whose length is bounded by the ranks
    int meld(int a, int b){
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        if (cmp(nodes[b].val, nodes[a].val)) swap(a, b);

        Node copy = nodes[a];
        nodes.push_back(copy);
        int c = (int)nodes.size() - 1;

        int right = meld(copy.right, b);  /// the call can reallocate nodes, assign after it
        nodes[c].right = right;
        if (rank(nodes[c].left) < rank(right)) swap(nodes[c].left, nodes[c].right);
        nodes[c].rank = rank(nodes[c].right) + 1;
        return c;
    }

    int rank(int h) const{
        return h == EMPTY ? 0 : nodes[h].rank;
    }
};
/// END SHARED LeftistHeap

const long long INF = LLONG_MAX;

struct Edge{
    int from, to, w;
};

/// Sidetrack: an edge off the shortest-path tree, keyed by the extra length it adds and remembering where it lands
using Sidetrack = pair<long long, int>;
using Heap = LeftistHeap<Sidetrack>;

int main(){
    int n, m, s, t, k;
    scanf("%d %d %d %d %d", &n, &m, &s, &t, &k);

    vector<Edge> edges(m);
    vector<vector<int>> out(n), in(n);
    for (int i = 0; i < m; i++){
        scanf("%d %d %d", &edges[i].from, &edges[i].to, &edges[i].w);
        out[edges[i].from].push_back(i);
        in[edges[i].to].push_back(i);
    }

    vector<long long> dis(n, INF);
    vector<int> par(n, -1), order;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;
    dis[t] = 0;
    pq.push({0, t});
    while (!pq.empty()){
        long long d = pq.top().first;
        int v = pq.top().second;
        pq.pop();
        if (d != dis[v]) continue;

        order.push_back(v);
        for (int i: in[v]){
            int u = edges[i].from;
            if (d + edges[i].w < dis[u]){
                dis[u] = d + edges[i].w;
                par[u] = i;
                pq.push({dis[u], u});
            }
        }
    }

    /// Dijkstra order puts every tree parent before its children, so each vertex's heap extends a finished one
    Heap heap;
    vector<int> root(n, Heap::EMPTY);
    for (int v: order){
        if (v != t) root[v] = root[edges[par[v]].to];
        for (int i: out[v]){
            int u = edges[i].to;
            if (i != par[v] && dis[u] != INF) root[v] = heap.push(root[v], {edges[i].w + dis[u] - dis[v], u});
        }
    }

    vector<long long> res;
    if (dis[s] != INF){
        res.push_back(dis[s]);
        if (root[s] != Heap::EMPTY) pq.push({dis[s] + heap.top(root[s]).first, root[s]});
    }

    /// A state is a heap node whose sidetrack was taken last: swap it for a child, or append the best sidetrack after it
    while ((int)res.size() < k && !pq.empty()){
        long long d = pq.top().first;
        int h = pq.top().second;
        pq.pop();
        res.push_back(d);

        long long without = d - heap.nodes[h].val.first;
        for (int c: {heap.nodes[h].left, heap.nodes[h].right}){
            if (c != Heap::EMPTY) pq.push({without + heap.nodes[c].val.first, c});
        }
        int next = root[heap.nodes[h].val.second];
        if (next != Heap::EMPTY) pq.push({d + heap.top(next).first, next});
    }

    res.resize(k, -1);
    for (long long x: res) printf("%lld\n", x);
    return 0;
}
