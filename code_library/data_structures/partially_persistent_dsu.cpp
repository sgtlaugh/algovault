/***
 *
 * Partially Persistent Disjoint Set Union
 * Unions at increasing times, connectivity and component size queries at any past time
 * Nodes can be 0 or 1 based, allocates nodes 0..n
 * connect returns true if a and b were in different components
 *
 * Complexity: O(log n) per operation, union by size without path compression, O(n) memory
 *
 * The k-th call to connect happens at time k, now is the number of calls so far
 * find_root(i, t), is_connected(a, b, t), component_size(i, t): the state after the first t calls, t >= 0
 * Any t >= now gives the current state
 * The first time a and b were connected is a binary search over t on is_connected
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PartiallyPersistentDSU{
    int now = 0;
    vector<int> parent, joined;  /// joined[i]: time i stopped being a root, INT_MAX while it still is one
    vector<vector<pair<int, int>>> sizes;  /// (time, size) of root i, times increasing

    PartiallyPersistentDSU(int n) : parent(n + 1), joined(n + 1, INT_MAX), sizes(n + 1, {{0, 1}}){
        iota(parent.begin(), parent.end(), 0);
    }

    int find_root(int i, int t){
        assert(t >= 0);
        while (i != parent[i] && joined[i] <= t) i = parent[i];
        return i;
    }

    bool connect(int a, int b){
        now++;
        a = find_root(a, now), b = find_root(b, now);
        if (a == b) return false;

        if (sizes[a].back().second > sizes[b].back().second) swap(a, b);
        parent[a] = b, joined[a] = now;
        sizes[b].push_back({now, sizes[a].back().second + sizes[b].back().second});
        return true;
    }

    bool is_connected(int a, int b, int t){
        return find_root(a, t) == find_root(b, t);
    }

    int component_size(int i, int t){
        auto& history = sizes[find_root(i, t)];
        return prev(upper_bound(history.begin(), history.end(), make_pair(t, INT_MAX)))->second;
    }
};

int main(){
    auto pp = PartiallyPersistentDSU(5);
    assert(pp.connect(1, 2) && pp.connect(3, 4) && !pp.connect(2, 1) && pp.connect(2, 4));
    assert(pp.now == 4);

    assert(!pp.is_connected(1, 2, 0) && pp.is_connected(1, 2, 1));
    assert(!pp.is_connected(1, 4, 3) && pp.is_connected(1, 4, 4) && pp.is_connected(3, 2, 100));
    assert(pp.find_root(1, 0) == 1 && pp.find_root(5, 4) == 5);

    assert(pp.component_size(1, 0) == 1 && pp.component_size(1, 1) == 2 && pp.component_size(4, 1) == 1);
    assert(pp.component_size(4, 3) == 2 && pp.component_size(3, 4) == 4 && pp.component_size(1, 100) == 4);
    assert(pp.component_size(5, 4) == 1);
    assert(pp.is_connected(3, 2, INT_MAX) && pp.component_size(1, INT_MAX) == 4 && !pp.is_connected(1, 5, INT_MAX));

    return 0;
}
