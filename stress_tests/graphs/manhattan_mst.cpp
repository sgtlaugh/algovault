#include "../common.h"

#define main library_main
#include "../../code_library/graphs/manhattan_mst.cpp"
#undef main

long long dist(const pair<long long, long long>& a, const pair<long long, long long>& b){
    return llabs(a.first - b.first) + llabs(a.second - b.second);
}

/// Kruskal on all n(n - 1) / 2 point pairs
long long brute_mst(const vector<pair<long long, long long>>& points){
    int n = points.size();
    vector<array<long long, 3>> edges;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++) edges.push_back({dist(points[i], points[j]), i, j});
    }
    sort(edges.begin(), edges.end());

    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };

    long long total = 0;
    for (auto [w, u, v] : edges){
        int a = find(u), b = find(v);
        if (a != b) parent[a] = b, total += w;
    }
    return total;
}

/// The result must be a spanning tree whose edge distances add up to the reported weight
void check(const vector<pair<long long, long long>>& points, long long expected){
    int n = points.size();
    vector<array<long long, 3>> candidates = manhattan_mst_candidates(points);
    assert((int)candidates.size() <= 4 * n);
    for (auto [w, u, v] : candidates) assert(u != v && w == dist(points[u], points[v]));

    ManhattanMSTResult r = manhattan_mst(points);
    assert(r.weight == expected);
    assert((int)r.edges.size() == max(0, n - 1));

    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    };

    long long sum = 0;
    for (auto [u, v] : r.edges){
        assert(0 <= u && u < n && 0 <= v && v < n);
        int a = find(u), b = find(v);
        assert(a != b);
        parent[a] = b, sum += dist(points[u], points[v]);
    }
    assert(sum == r.weight);
}

int main(){
    const long long LIM = 1000000000;

    /// Every multiset of up to 4 points on a 3 x 3 grid, as ordered sequences so input order varies too
    for (int n = 0; n <= 4; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 9;
        for (int mask = 0; mask < total; mask++){
            vector<pair<long long, long long>> points;
            for (int i = 0, m = mask; i < n; i++, m /= 9) points.push_back({m % 9 / 3, m % 3});
            check(points, brute_mst(points));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 300 ? it % 30 + 1 : stress::rand_int(1, 120);
        long long range = it % 4 == 0 ? 3 : it % 4 == 1 ? 20 : it % 4 == 2 ? 1000 : LIM;
        vector<pair<long long, long long>> points(n);
        for (auto& [x, y] : points){
            x = stress::rand_int(-range, range), y = stress::rand_int(-range, range);
            if (it % 7 == 0) x = stress::rand_int(0, 1) ? -range : range;
            if (it % 11 == 0) y = max(-range, min(range, x + stress::rand_int(-1, 1)));
        }
        check(points, brute_mst(points));
    }

    /// A large shuffled grid with spacing s: every point has a neighbour at distance s, so the MST weighs (n - 1) * s
    int side = 320, n = side * side;
    long long s = 2 * LIM / (side - 1);
    vector<pair<long long, long long>> points;
    for (int i = 0; i < side; i++){
        for (int j = 0; j < side; j++) points.push_back({-LIM + i * s, -LIM + j * s});
    }
    shuffle(points.begin(), points.end(), stress::rng());
    check(points, (n - 1) * s);
    return 0;
}
