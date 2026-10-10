// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/k_shortest_walk
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/leftist_heap.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    const long long INF = LLONG_MAX;
    int n = read_int(), m = read_int(), s = read_int(), t = read_int(), k = read_int();
    vector<int> from(m), to(m);
    vector<long long> cost(m);
    vector<vector<int>> out(n), in(n);
    for (int i = 0; i < m; i++){
        from[i] = read_int(), to[i] = read_int(), cost[i] = read_int();
        out[from[i]].push_back(i);
        in[to[i]].push_back(i);
    }

    vector<long long> dist(n, INF);
    vector<int> tree_edge(n, -1), order;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
    dist[t] = 0;
    pq.push({0, t});
    while (!pq.empty()){
        auto [d, v] = pq.top();
        pq.pop();
        if (d != dist[v]) continue;

        order.push_back(v);
        for (int e: in[v]){
            if (d + cost[e] < dist[from[e]]){
                dist[from[e]] = d + cost[e], tree_edge[from[e]] = e;
                pq.push({dist[from[e]], from[e]});
            }
        }
    }

    /// Eppstein: a walk is its sequence of sidetracks, edges off the shortest path tree towards t
    /// Vertex v's heap holds (detour cost, head) of the sidetracks leaving v and its tree ancestors
    using Heap = LeftistHeap<pair<long long, int>>;
    Heap heap;
    vector<int> root(n, Heap::EMPTY);
    for (int v: order){
        int h = (v == t) ? Heap::EMPTY : root[to[tree_edge[v]]];
        for (int e: out[v]){
            if (e != tree_edge[v] && dist[to[e]] != INF) h = heap.push(h, {cost[e] + dist[to[e]] - dist[v], to[e]});
        }
        root[v] = h;
    }

    vector<long long> res;
    if (dist[s] != INF){
        res.push_back(dist[s]);

        /// (walk length, heap node): the walk takes the sidetracks leading to the node, then the node's own
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> walks;
        if (root[s] != Heap::EMPTY) walks.push({dist[s] + heap.top(root[s]).first, root[s]});
        while ((int)res.size() < k && !walks.empty()){
            auto [d, h] = walks.top();
            walks.pop();
            res.push_back(d);

            auto node = heap.nodes[h];
            for (int child: {node.left, node.right}){
                if (child != Heap::EMPTY) walks.push({d - node.val.first + heap.top(child).first, child});
            }

            int next = root[node.val.second];
            if (next != Heap::EMPTY) walks.push({d + heap.top(next).first, next});
        }
    }

    for (int i = 0; i < k; i++) printf("%lld\n", i < (int)res.size() ? res[i] : -1LL);
    return 0;
}
