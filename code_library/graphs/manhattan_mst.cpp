/***
 *
 * Manhattan Minimum Spanning Tree
 * Minimum spanning tree of points in the plane under the L1 distance |x1 - x2| + |y1 - y2|
 *
 * Complexity: O(n log n)
 *
 * manhattan_mst_candidates(points): at most 4n edges {distance, i, j} that contain a minimum spanning tree
 *     one sweep per octant pair keeps, for every point, only its nearest neighbour in that octant
 *     feed them to any MST algorithm when the points are part of a larger graph
 * manhattan_mst(points): runs Kruskal on the candidates, returns {total weight, chosen edges {i, j}}
 *     n - 1 edges for n >= 1 points, none for an empty input; duplicate points are joined by zero-weight edges
 *
 * points[i] = {x, y} with |x|, |y| <= 1e9, so every distance and the total weight of up to 1e9 edges fit in long long
 *
 * Example:
 *     ManhattanMSTResult r = manhattan_mst({{0, 0}, {2, 2}, {5, 2}});
 *     // r.weight = 7, r.edges joins points 1-2 and 0-1 (each pair in either orientation)
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
    /***
     *   A(0, 0), B(2, 2), C(3, 10), D(5, 2), E(7, 0)
     *   BD = 3, AB = 4, DE = 4, BC = 9 beat AD = AE = BE = 7 once A, B, D, E are joined
    ***/
    ManhattanMSTResult r = manhattan_mst({{0, 0}, {2, 2}, {3, 10}, {5, 2}, {7, 0}});
    assert(r.weight == 20 && r.edges.size() == 4);
    set<set<int>> chosen;
    for (auto [u, v] : r.edges) chosen.insert({u, v});
    assert((chosen == set<set<int>>{{1, 3}, {0, 1}, {3, 4}, {1, 2}}));

    assert(manhattan_mst({{0, 0}, {2, 2}, {5, 2}}).weight == 7);
    assert(manhattan_mst({{1, 1}, {1, 1}, {4, 5}}).weight == 7);
    assert(manhattan_mst({{0, 0}, {0, 5}, {0, 2}}).weight == 5);
    assert(manhattan_mst({{3, -1}, {-2, -1}, {7, -1}, {-2, -1}}).weight == 9);
    assert(manhattan_mst({{-1000000000, -1000000000}, {1000000000, 1000000000}}).weight == 4000000000LL);
    assert(manhattan_mst({{1000000000, -1000000000}, {-1000000000, 1000000000}, {0, 0}}).weight == 4000000000LL);

    assert(manhattan_mst({{7, 7}}).weight == 0 && manhattan_mst({{7, 7}}).edges.empty());
    assert(manhattan_mst({}).weight == 0 && manhattan_mst({}).edges.empty());
    assert(manhattan_mst_candidates({}).empty());
    return 0;
}
