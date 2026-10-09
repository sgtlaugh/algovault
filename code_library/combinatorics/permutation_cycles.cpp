/***
 *
 * Permutation Cycles
 * Cycle decomposition, k-th power, k-th root, parity and minimum sorting swaps of a permutation
 *
 * Complexity: O(n) per call, plus O(sqrt(n) log^2 k) for permutation_root
 *
 * p is a permutation of 0..n-1, p^k applies p k times: p^2[i] = p[p[i]]
 *
 * permutation_cycles(p): the cycles, each listed as i, p[i], p[p[i]], ...
 * permutation_power(p, k): p^k for any long long k, negative k gives powers of the inverse
 * permutation_root(p, k): some q with q^k = p for 1 <= k <= LLONG_MAX, or {} when none exists (n >= 1)
 *     among all k-th roots it returns one with the most cycles, so the fewest sorting swaps
 * sorting_swaps(p): the minimum list of position swaps that turns the array p into the identity, n - cycles long
 *     replaying it in reverse on the identity rebuilds p
 * permutation_parity(p): 0 for even, 1 for odd, (n - cycles) mod 2
 *
 * Root: q^k splits a q-cycle of length m into gcd(m, k) cycles of length m / gcd(m, k)
 * So the cycles of p of length len are merged in groups of t, the smallest t with gcd(t * len, k) = t
 * No root exists when the number of len-cycles is not a multiple of t
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<vector<int>> permutation_cycles(const vector<int>& p){
    int n = p.size();
    vector<vector<int>> cycles;
    vector<bool> visited(n, false);

    for (int i = 0; i < n; i++){
        if (visited[i]) continue;
        cycles.emplace_back();
        for (int j = i; !visited[j]; j = p[j]){
            visited[j] = true;
            cycles.back().push_back(j);
        }
    }
    return cycles;
}

int permutation_parity(const vector<int>& p){
    return ((int)p.size() - (int)permutation_cycles(p).size()) & 1;
}

vector<int> permutation_power(const vector<int>& p, long long k){
    vector<int> res(p.size());
    for (const auto& cycle: permutation_cycles(p)){
        int len = cycle.size();
        int shift = (k % len + len) % len;
        for (int j = 0; j < len; j++) res[cycle[j]] = cycle[(j + shift) % len];
    }
    return res;
}

vector<int> permutation_root(const vector<int>& p, long long k){
    assert(k >= 1);
    int n = p.size();
    vector<int> res(n);
    vector<vector<vector<int>>> by_len(n + 1);
    for (auto& cycle: permutation_cycles(p)) by_len[cycle.size()].push_back(move(cycle));

    for (int len = 1; len <= n; len++){
        const auto& cycles = by_len[len];
        if (cycles.empty()) continue;

        long long t = 1, rest = k;
        for (long long g = __gcd((long long)len, rest); g != 1; g = __gcd((long long)len, rest)){
            t *= g;
            rest /= g;
        }
        if ((long long)cycles.size() % t != 0) return {};

        int m = t * len;
        long long step = k % m;
        vector<int> merged(m);
        for (int i = 0; i < (int)cycles.size(); i += t){
            for (int x = 0; x < t; x++){
                for (int y = 0; y < len; y++) merged[(x + y * step) % m] = cycles[i + x][y];
            }
            for (int j = 0; j < m; j++) res[merged[j]] = merged[(j + 1) % m];
        }
    }
    return res;
}

vector<pair<int, int>> sorting_swaps(vector<int> p){
    vector<pair<int, int>> swaps;
    for (int i = 0; i < (int)p.size(); i++){
        while (p[i] != i){
            swaps.emplace_back(i, p[i]);
            swap(p[i], p[p[i]]);
        }
    }
    return swaps;
}

int main(){
    vector<int> p = {1, 2, 0, 4, 3};
    assert((permutation_cycles(p) == vector<vector<int>>{{0, 1, 2}, {3, 4}}));
    assert((permutation_power(p, 0) == vector<int>{0, 1, 2, 3, 4}));
    assert((permutation_power(p, 2) == vector<int>{2, 0, 1, 3, 4}));
    assert((permutation_power(p, 6) == vector<int>{0, 1, 2, 3, 4}));
    assert((permutation_power(p, -1) == vector<int>{2, 0, 1, 4, 3}));
    assert((permutation_power(p, 1000000000000000000LL) == vector<int>{1, 2, 0, 3, 4}));
    assert((permutation_power(p, LLONG_MIN) == vector<int>{1, 2, 0, 3, 4}));
    assert(permutation_parity(p) == 1);
    assert((sorting_swaps(p) == vector<pair<int, int>>{{0, 1}, {0, 2}, {3, 4}}));

    vector<int> two_swaps = {1, 0, 3, 2};
    assert((permutation_root(two_swaps, 2) == vector<int>{2, 3, 1, 0}));
    assert(permutation_parity(two_swaps) == 0);
    assert((permutation_root({1, 0}, 2) == vector<int>{}));
    assert((permutation_root({1, 2, 0}, 2) == vector<int>{2, 0, 1}));
    assert((permutation_root({1, 2, 0}, 3) == vector<int>{}));
    assert((permutation_root({1, 2, 0}, 4) == vector<int>{1, 2, 0}));
    assert((permutation_root({0, 1, 2}, LLONG_MAX) == vector<int>{0, 1, 2}));
    assert((permutation_root({1, 0}, 1LL << 62) == vector<int>{}));
    assert((permutation_root({1, 0}, 3) == vector<int>{1, 0}));

    assert(permutation_cycles({}).empty() && permutation_power({}, 5).empty() && sorting_swaps({}).empty());
    assert((permutation_power({0}, -7) == vector<int>{0}) && permutation_parity({0}) == 0);

    return 0;
}
