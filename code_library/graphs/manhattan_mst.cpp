/***
 *
 * Manhattan Minimum Spanning Tree
 * Minimum spanning tree of points in the plane under the L1 distance |x1 - x2| + |y1 - y2|
 *
 * Complexity: O(n log n)
 *
 * manhattan_mst(points): {weight, edges}, the total weight and the n - 1 chosen edges {i, j}, none for an empty input
 * manhattan_mst_candidates(points): at most 4n edges {distance, i, j} that contain a minimum spanning tree
 *
 * Notes:
 *   - points[i] = {x, y} with |x|, |y| <= 1e9, so every distance and the total weight fit in long long
 *   - Duplicate points are joined by zero-weight edges
 *   - Feed the candidates to any MST algorithm when the points are part of a larger graph
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ManhattanMSTResult{
    long long weight;
    vector<array<int, 2>> edges;
};

vector<array<long long, 3>> manhattan_mst_candidates(vector<pair<long long, long long>> points){
    int n = points.size();
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    vector<array<long long, 3>> edges;
    edges.reserve(4 * n);

    for (int k = 0; k < 4; k++){
        sort(order.begin(), order.end(), [&](int i, int j){
            return points[i].first + points[i].second < points[j].first + points[j].second;
        });

        /// Keyed by -y: the points still waiting for their nearest neighbour in this octant
        map<long long, int> sweep;
        for (int i : order){
            for (auto it = sweep.lower_bound(-points[i].second); it != sweep.end(); sweep.erase(it++)){
                int j = it->second;
                long long dx = points[i].first - points[j].first, dy = points[i].second - points[j].second;
                if (dy > dx) break;
                edges.push_back({dx + dy, i, j});
            }
            sweep[-points[i].second] = i;
        }

        /// Reflect and rotate so the next sweep covers the next octant
        for (auto& p : points){
            if (k & 1) p.first = -p.first;
            else swap(p.first, p.second);
        }
    }
    return edges;
}

ManhattanMSTResult manhattan_mst(const vector<pair<long long, long long>>& points){
    int n = points.size();
    vector<array<long long, 3>> edges = manhattan_mst_candidates(points);
    sort(edges.begin(), edges.end());

    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    };

    ManhattanMSTResult res{0, {}};
    for (auto [w, u, v] : edges){
        int a = find(u), b = find(v);
        if (a == b) continue;
        parent[a] = b;
        res.weight += w;
        res.edges.push_back({(int)u, (int)v});
    }
    return res;
}

int main(){
    /// Naming the points A to E in order, AB = 4, BD = 3, DE = 4 and BC = 9 beat AD = AE = BE = 7
    ManhattanMSTResult r = manhattan_mst({{0, 0}, {2, 2}, {3, 10}, {5, 2}, {7, 0}});
    assert(r.weight == 20);
    set<set<int>> chosen;
    for (auto [u, v] : r.edges) chosen.insert({u, v});
    assert((chosen == set<set<int>>{{0, 1}, {1, 3}, {3, 4}, {1, 2}}));
    return 0;
}
