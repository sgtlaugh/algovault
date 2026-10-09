#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/permutation_tree.cpp"
#undef main

long long brute_count(const vector<int>& p, vector<vector<char>>* good = nullptr){
    int n = p.size();
    long long total = 0;
    for (int l = 0; l < n; l++){
        int lo = p[l], hi = p[l];
        for (int r = l; r < n; r++){
            lo = min(lo, p[r]), hi = max(hi, p[r]);
            if (hi - lo != r - l) continue;
            total++;
            if (good) (*good)[l][r] = 1;
        }
    }
    return total;
}

void check_tree(const vector<int>& p){
    int n = p.size();
    PermutationTree pt(p);
    vector<vector<char>> good(n, vector<char>(n)), from_tree(n, vector<char>(n));
    assert(pt.count_intervals() == brute_count(p, &good));
    if (n == 0){
        assert(pt.root == -1);
        return;
    }

    int nodes = pt.children.size();
    assert(nodes <= 2 * n - 1);
    assert(pt.parent[pt.root] == -1 && pt.span[pt.root] == make_pair(0, n - 1));
    for (int v = 0; v < nodes; v++){
        auto [l, r] = pt.span[v];
        int lo = *min_element(p.begin() + l, p.begin() + r + 1), hi = *max_element(p.begin() + l, p.begin() + r + 1);
        assert(pt.range[v] == make_pair(lo, hi) && hi - lo == r - l);
        assert(v == pt.root || pt.parent[v] != -1);

        const auto& ch = pt.children[v];
        int k = ch.size();
        if (v < n){
            assert(k == 0 && l == v && r == v && pt.type[v] == PermutationTree::CUT);
            from_tree[l][r] = 1;
            continue;
        }

        assert(pt.span[ch[0]].first == l && pt.span[ch[k - 1]].second == r);
        for (int i = 0; i < k; i++) assert(pt.parent[ch[i]] == v);
        for (int i = 0; i + 1 < k; i++) assert(pt.span[ch[i]].second + 1 == pt.span[ch[i + 1]].first);

        if (pt.type[v] == PermutationTree::CUT){
            assert(k >= 4);
            from_tree[l][r] = 1;
            continue;
        }

        assert(k >= 2);
        for (int i = 0; i < k; i++) assert(pt.type[ch[i]] != pt.type[v]);
        for (int i = 0; i + 1 < k; i++){
            if (pt.type[v] == PermutationTree::INCREASING) assert(pt.range[ch[i]].second + 1 == pt.range[ch[i + 1]].first);
            else assert(pt.range[ch[i + 1]].second + 1 == pt.range[ch[i]].first);
        }
        for (int i = 0; i < k; i++){
            for (int j = i + 1; j < k; j++) from_tree[pt.span[ch[i]].first][pt.span[ch[j]].second] = 1;
        }
    }
    assert(from_tree == good);
}

vector<int> random_permutation(int n){
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    shuffle(p.begin(), p.end(), stress::rng());
    return p;
}

/// Shuffles small blocks of a sorted or reversed array so join and cut nodes nest several levels deep
vector<int> blocky_permutation(int n){
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    if (stress::rand_int(0, 1)) reverse(p.begin(), p.end());
    for (int rounds = stress::rand_int(0, 4); rounds > 0; rounds--){
        int len = stress::rand_int(1, max(1, n / 2)), start = stress::rand_int(0, n - len);
        shuffle(p.begin() + start, p.begin() + start + len, stress::rng());
    }
    return p;
}

/// O(n^2) scan marks every [l, r] with max - min == r - l, then checks the tree's node invariants and that
/// leaves, whole cut nodes and runs of 2+ join children are exactly those intervals
int main(){
    for (int n = 0; n <= 7; n++){
        vector<int> p(n);
        iota(p.begin(), p.end(), 0);
        do check_tree(p); while (next_permutation(p.begin(), p.end()));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it < 2000 ? 40 : 150);
        check_tree(it % 2 ? random_permutation(n) : blocky_permutation(n));
    }

    for (int n : {3000, 3001}){
        for (int kind = 0; kind < 2; kind++){
            auto p = kind ? random_permutation(n) : blocky_permutation(n);
            assert(PermutationTree(p).count_intervals() == brute_count(p));
        }
    }

    int big = 200000;
    vector<int> up(big);
    iota(up.begin(), up.end(), 0);
    assert(PermutationTree(up).count_intervals() == 1LL * big * (big + 1) / 2);
    reverse(up.begin(), up.end());
    assert(PermutationTree(up).count_intervals() == 1LL * big * (big + 1) / 2);
    PermutationTree shuffled(random_permutation(big));
    assert((int)shuffled.children.size() <= 2 * big - 1 && shuffled.span[shuffled.root] == make_pair(0, big - 1));

    return 0;
}
