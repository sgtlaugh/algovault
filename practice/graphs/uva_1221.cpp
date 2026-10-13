/***
 *
 * Hopcroft Karp
 * Maximum matching in a bipartite graph
 *
 * Complexity: O(m * sqrt(n)), O(n + m) memory
 *
 * HopcroftKarp hk(n_left, n_right): left vertices 0..n_left - 1, right vertices 0..n_right - 1
 * hk.add_edge(u, v): edge from left u to right v
 * hk.max_matching(): size of a maximum matching, can be called again after adding more edges
 * hk.match_left[u]: right vertex matched to u or -1, hk.match_right[v]: left vertex matched to v or -1
 *
 * Each phase augments along shortest paths only, which is what bounds the phases by O(sqrt(n))
 * The augmenting search is iterative, so long augmenting paths cannot overflow the stack
 *
 * minimum_path_cover.cpp carries a checked copy of this struct
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED hopcroft_karp
struct HopcroftKarp{
    int n_left, n_right;
    vector<vector<int>> adj;
    vector<int> match_left, match_right, dist, it;

    HopcroftKarp(int n_left, int n_right) : n_left(n_left), n_right(n_right), adj(n_left), match_left(n_left, -1), match_right(n_right, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    /// Layers left vertices by distance from free left vertices, stopping at the first layer that reaches a free right vertex
    int bfs(){
        dist.assign(n_left, -1);
        vector<int> queue;
        for (int u = 0; u < n_left; u++){
            if (match_left[u] == -1) dist[u] = 0, queue.push_back(u);
        }

        int limit = INT_MAX;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int v : adj[u]){
                int w = match_right[v];
                if (w == -1) limit = min(limit, dist[u]);
                else if (dist[w] == -1 && dist[u] < limit){
                    dist[w] = dist[u] + 1;
                    queue.push_back(w);
                }
            }
        }
        return limit;
    }

    void augment(int root, int limit){
        vector<int> stack = {root};
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == (int)adj[u].size()){
                dist[u] = -1;
                stack.pop_back();
                continue;
            }

            int v = adj[u][it[u]++], w = match_right[v];
            if (w == -1 && dist[u] == limit){
                for (int x : stack){
                    int y = adj[x][it[x] - 1];
                    match_left[x] = y, match_right[y] = x;
                }
                return;
            }
            if (w != -1 && dist[w] == dist[u] + 1) stack.push_back(w);
        }
    }

    int max_matching(){
        while (true){
            int limit = bfs();
            if (limit == INT_MAX) break;
            it.assign(n_left, 0);
            for (int u = 0; u < n_left; u++){
                if (match_left[u] == -1) augment(u, limit);
            }
        }
        return n_left - count(match_left.begin(), match_left.end(), -1);
    }
};
/// END SHARED hopcroft_karp

int main(){
    int h, a;
    while (scanf("%d %d", &h, &a) == 2 && (h || a)){
        vector<long long> hn(h), hp(h), an(a), ap(a);
        for (int i = 0; i < h; i++) scanf("%lld %lld", &hn[i], &hp[i]);
        for (int j = 0; j < a; j++) scanf("%lld %lld", &an[j], &ap[j]);

        /// win[i][j]: earliest arrival year at which planet i's fleet beats alien j, leaving at year arrival - travel >= 0
        const long long NEVER = LLONG_MAX;
        vector<vector<long long>> win(h, vector<long long>(a));
        vector<long long> years;
        for (int i = 0; i < h; i++){
            for (int j = 0; j < a; j++){
                long long t, gap;
                scanf("%lld", &t);
                gap = an[j] + t * ap[j] - hn[i];
                if (gap <= 0) win[i][j] = t;
                else if (hp[i] > ap[j]) win[i][j] = t + (gap + hp[i] - ap[j] - 1) / (hp[i] - ap[j]);
                else win[i][j] = NEVER;
                if (win[i][j] != NEVER) years.push_back(win[i][j]);
            }
        }
        sort(years.begin(), years.end());
        years.erase(unique(years.begin(), years.end()), years.end());

        /// Smallest year by which every alien planet can get its own human planet
        auto feasible = [&](long long year){
            HopcroftKarp hk(a, h);
            for (int i = 0; i < h; i++){
                for (int j = 0; j < a; j++){
                    if (win[i][j] <= year) hk.add_edge(j, i);
                }
            }
            return hk.max_matching() == a;
        };

        int lo = 0, hi = years.size();
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (feasible(years[mid])) hi = mid;
            else lo = mid + 1;
        }
        if (lo == (int)years.size()) puts("IMPOSSIBLE");
        else printf("%lld\n", years[lo]);
    }
    return 0;
}
