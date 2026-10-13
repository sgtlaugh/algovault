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
 * find(i, t), same(a, b, t), size(i, t): the state after the first t calls, t >= 0
 * Any t >= now gives the current state
 * The first time a and b were connected is a binary search over t on same
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

    int find(int i, int t){
        assert(t >= 0);
        while (i != parent[i] && joined[i] <= t) i = parent[i];
        return i;
    }

    bool connect(int a, int b){
        now++;
        a = find(a, now), b = find(b, now);
        if (a == b) return false;

        if (sizes[a].back().second > sizes[b].back().second) swap(a, b);
        parent[a] = b, joined[a] = now;
        sizes[b].push_back({now, sizes[a].back().second + sizes[b].back().second});
        return true;
    }

    bool same(int a, int b, int t){
        return find(a, t) == find(b, t);
    }

    int size(int i, int t){
        auto& history = sizes[find(i, t)];
        return prev(upper_bound(history.begin(), history.end(), make_pair(t, INT_MAX)))->second;
    }
};

int main(){
    auto pp = PartiallyPersistentDSU(5);
    assert(pp.connect(1, 2) && pp.connect(3, 4) && !pp.connect(2, 1) && pp.connect(2, 4));
    assert(pp.now == 4);

    assert(!pp.same(1, 2, 0) && pp.same(1, 2, 1));
    assert(!pp.same(1, 4, 3) && pp.same(1, 4, 4) && pp.same(3, 2, 100));
    assert(pp.find(1, 0) == 1 && pp.find(5, 4) == 5);

    assert(pp.size(1, 0) == 1 && pp.size(1, 1) == 2 && pp.size(4, 1) == 1);
    assert(pp.size(4, 3) == 2 && pp.size(3, 4) == 4 && pp.size(1, 100) == 4);
    assert(pp.size(5, 4) == 1);
    assert(pp.same(3, 2, INT_MAX) && pp.size(1, INT_MAX) == 4 && !pp.same(1, 5, INT_MAX));

    return 0;
}
