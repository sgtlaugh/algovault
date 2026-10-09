#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/permutation_cycles.cpp"
#undef main

vector<int> compose(const vector<int>& a, const vector<int>& b){
    vector<int> res(a.size());
    for (int i = 0; i < (int)a.size(); i++) res[i] = a[b[i]];
    return res;
}

/// p^k by binary exponentiation of composition, negative k through the inverse
vector<int> brute_power(vector<int> p, long long k){
    int n = p.size();
    vector<int> res(n);
    iota(res.begin(), res.end(), 0);
    if (k < 0){
        vector<int> inverse(n);
        for (int i = 0; i < n; i++) inverse[p[i]] = i;
        p = inverse;
    }

    for (unsigned long long e = k < 0 ? -(unsigned long long)k : k; e; e >>= 1){
        if (e & 1) res = compose(res, p);
        p = compose(p, p);
    }
    return res;
}

/// A cycle is counted at its smallest element
int brute_cycle_count(const vector<int>& p){
    int count = 0;
    for (int i = 0; i < (int)p.size(); i++){
        int j = p[i];
        while (j > i) j = p[j];
        count += j == i;
    }
    return count;
}

int inversion_parity(const vector<int>& p){
    int parity = 0;
    for (int i = 0; i < (int)p.size(); i++){
        for (int j = i + 1; j < (int)p.size(); j++) parity ^= p[i] > p[j];
    }
    return parity;
}

void check_cycles(const vector<int>& p){
    int n = p.size();
    vector<int> seen(n, 0);
    auto cycles = permutation_cycles(p);
    assert((int)cycles.size() == brute_cycle_count(p));
    for (const auto& cycle: cycles){
        int len = cycle.size();
        for (int j = 0; j < len; j++){
            seen[cycle[j]]++;
            assert(p[cycle[j]] == cycle[(j + 1) % len]);
        }
    }
    for (int i = 0; i < n; i++) assert(seen[i] == 1);
}

void check_swaps(const vector<int>& p, int expected_count){
    auto swaps = sorting_swaps(p);
    assert((int)swaps.size() == expected_count);
    vector<int> a = p;
    for (auto [x, y]: swaps) swap(a[x], a[y]);
    for (int i = 0; i < (int)a.size(); i++) assert(a[i] == i);

    vector<int> b(p.size());
    iota(b.begin(), b.end(), 0);
    for (int i = (int)swaps.size() - 1; i >= 0; i--) swap(b[swaps[i].first], b[swaps[i].second]);
    assert(b == p);
}

void check_root(const vector<int>& p, long long k, bool exists){
    auto root = permutation_root(p, k);
    if (!exists){
        assert(root.empty());
        return;
    }
    assert(root.size() == p.size() && brute_power(root, k) == p);
}

/// Most cycles of any k-th root of p, or -1 when none exists
/// A root cycle of length s * len splits under q^k into len-cycles exactly when gcd(s * len, k) = s,
/// so per length a knapsack partitions the count of len-cycles into such group sizes s, maximizing groups
int root_cycles_oracle(const vector<int>& p, long long k){
    int n = p.size();
    vector<int> count(n + 1, 0);
    for (const auto& cycle: permutation_cycles(p)) count[cycle.size()]++;

    int total = 0;
    for (int len = 1; len <= n; len++){
        int c = count[len];
        if (c == 0) continue;
        vector<int> sizes;
        for (int s = 1; s <= c; s++){
            if (gcd((long long)s * len, k) == s) sizes.push_back(s);
        }

        vector<int> best(c + 1, -1);
        best[0] = 0;
        for (int x = 1; x <= c; x++){
            for (int s: sizes){
                if (s <= x && best[x - s] >= 0) best[x] = max(best[x], best[x - s] + 1);
            }
        }
        if (best[c] < 0) return -1;
        total += best[c];
    }
    return total;
}

void check_root_against_oracle(const vector<int>& p, long long k){
    int expected = root_cycles_oracle(p, k);
    check_root(p, k, expected >= 0);
    if (expected >= 0) assert(brute_cycle_count(permutation_root(p, k)) == expected);
}

/// Every permutation of n <= 7 against all k-th roots found by powering every permutation
void exhaustive(int n){
    vector<vector<int>> perms;
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    do perms.push_back(p);
    while (next_permutation(p.begin(), p.end()));

    /// Minimum swaps by BFS from the identity, a swap is its own inverse so distance(identity, p) = distance(p, identity)
    map<vector<int>, int> dist = {{perms[0], 0}};
    queue<vector<int>> q;
    q.push(perms[0]);
    while (!q.empty()){
        auto cur = q.front();
        q.pop();
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                auto next = cur;
                swap(next[i], next[j]);
                if (dist.count(next)) continue;
                dist[next] = dist[cur] + 1;
                q.push(next);
            }
        }
    }

    for (const auto& perm: perms){
        check_cycles(perm);
        check_swaps(perm, dist[perm]);
        assert(permutation_parity(perm) == inversion_parity(perm));
        for (long long k = -13; k <= 13; k++) assert(permutation_power(perm, k) == brute_power(perm, k));
    }

    for (long long k = 1; k <= 13; k++){
        map<vector<int>, int> most_cycles;
        for (const auto& root: perms){
            auto target = brute_power(root, k);
            auto it = most_cycles.find(target);
            int cycles = brute_cycle_count(root);
            if (it == most_cycles.end()) most_cycles[target] = cycles;
            else it->second = max(it->second, cycles);
        }
        for (const auto& perm: perms){
            auto it = most_cycles.find(perm);
            check_root(perm, k, it != most_cycles.end());
            if (it != most_cycles.end()) assert(brute_cycle_count(permutation_root(perm, k)) == it->second);
            assert(root_cycles_oracle(perm, k) == (it == most_cycles.end() ? -1 : it->second));
        }
    }
}

vector<int> random_permutation(int n){
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    shuffle(p.begin(), p.end(), stress::rng());
    return p;
}

/// Many short cycles so that k-th powers split and roots must merge several cycles of one length
vector<int> short_cycles_permutation(int n, int max_len){
    vector<int> order = random_permutation(n), p(n);
    for (int i = 0; i < n;){
        int len = min<int>(n - i, stress::rand_int(1, max_len));
        for (int j = 0; j < len; j++) p[order[i + j]] = order[i + (j + 1) % len];
        i += len;
    }
    return p;
}

long long random_exponent(){
    const long long special[] = {1, 2, 4, 6, 12, 60, 720, 1LL << 62, 1000000000000000000LL, LLONG_MAX, 2LL * 3 * 5 * 7 * 11 * 13};
    int kind = stress::rand_int(0, 2);
    if (kind == 0) return special[stress::rand_int(0, 10)];
    if (kind == 1) return stress::rand_int(1, 64);
    return stress::rand_int(1, LLONG_MAX);
}

int main(){
    for (int n = 0; n <= 7; n++) exhaustive(n);

    for (long long it = 0; it < stress::scaled(400); it++){
        int n = it < 100 ? it + 1 : stress::rand_int(1, 3000);
        vector<int> p = it % 2 ? random_permutation(n) : short_cycles_permutation(n, stress::rand_int(1, 12));

        check_cycles(p);
        check_swaps(p, n - brute_cycle_count(p));
        if (n <= 400) assert(permutation_parity(p) == inversion_parity(p));
        assert(permutation_parity(p) == (n - brute_cycle_count(p)) % 2);

        long long k = random_exponent();
        if (stress::rand_int(0, 3) == 0) k = stress::rand_int(0, 1) ? -k : LLONG_MIN;
        assert(permutation_power(p, k) == brute_power(p, k));

        k = random_exponent();
        auto planted = brute_power(p, k);
        assert(root_cycles_oracle(planted, k) >= 0);
        check_root_against_oracle(planted, k);
        check_root_against_oracle(p, k);
    }

    return 0;
}
