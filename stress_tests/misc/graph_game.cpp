#include "../common.h"

#define main library_main
#include "../../code_library/misc/graph_game.cpp"
#undef main

const int BIG = 1 << 20;

/// Minimax over a horizon of 2n + 2 moves, a value BIG - d wins in d moves, d - BIG loses in d moves, 0 draws
vector<int> brute(const vector<vector<int>>& adj, bool no_move_wins){
    int n = adj.size();
    int stuck = no_move_wins ? BIG : -BIG;
    vector<int> val(n);
    for (int u = 0; u < n; u++) val[u] = adj[u].empty() ? stuck : 0;

    for (int round = 0; round < 2 * n + 2; round++){
        vector<int> nxt(n);
        for (int u = 0; u < n; u++){
            if (adj[u].empty()){
                nxt[u] = stuck;
                continue;
            }
            nxt[u] = INT_MIN;
            for (int v : adj[u]){
                int mine = val[v] > 0 ? 1 - val[v] : val[v] < 0 ? -val[v] - 1 : 0;
                nxt[u] = max(nxt[u], mine);
            }
        }
        val = nxt;
    }

    return val;
}

void check(const vector<vector<int>>& adj){
    int n = adj.size();
    GraphGame g(adj);
    vector<int> val = brute(adj, false);
    for (int u = 0; u < n; u++){
        if (val[u] > 0) assert(g.result[u] == GraphGame::WIN && g.moves[u] == BIG - val[u]);
        if (val[u] < 0) assert(g.result[u] == GraphGame::LOSE && g.moves[u] == val[u] + BIG);
        if (val[u] == 0) assert(g.result[u] == GraphGame::DRAW && g.moves[u] == -1);
    }

    /// The header's recipe for "a player who cannot move wins": one extra state, every finite count 1 higher
    vector<vector<int>> with_sink = adj;
    with_sink.push_back({});
    for (int u = 0; u < n; u++) if (adj[u].empty()) with_sink[u].push_back(n);
    GraphGame flipped(with_sink);
    vector<int> flipped_val = brute(adj, true);
    for (int u = 0; u < n; u++){
        if (flipped_val[u] > 0) assert(flipped.result[u] == GraphGame::WIN && flipped.moves[u] == BIG - flipped_val[u] + 1);
        if (flipped_val[u] < 0) assert(flipped.result[u] == GraphGame::LOSE && flipped.moves[u] == flipped_val[u] + BIG + 1);
        if (flipped_val[u] == 0) assert(flipped.result[u] == GraphGame::DRAW);
    }
}

vector<vector<int>> random_game(int n, int shape){
    vector<vector<int>> adj(n);
    if (shape == 0){
        int m = stress::rand_int(0, 3 * n);
        for (int i = 0; i < m; i++) adj[stress::rand_int(0, n - 1)].push_back(stress::rand_int(0, n - 1));
    }
    else{
        /// Mostly backward moves with a few forward ones, so long forced lines meet a few cycles
        for (int u = 1; u < n; u++){
            int k = stress::rand_int(1, 3);
            for (int i = 0; i < k; i++) adj[u].push_back(stress::rand_int(max(0, u - 4), u - 1));
            if (stress::rand_int(0, 9) == 0) adj[u].push_back(stress::rand_int(0, n - 1));
        }
    }

    return adj;
}

int main(){
    GraphGame empty(vector<vector<int>>{});
    assert(empty.result.empty());

    /// Every directed graph on up to 4 states, self loops included
    for (int n = 1; n <= 4; n++){
        for (int mask = 0; mask < (1 << (n * n)); mask++){
            vector<vector<int>> adj(n);
            for (int b = 0; b < n * n; b++) if (mask >> b & 1) adj[b / n].push_back(b % n);
            check(adj);
        }
    }

    for (int it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 40);
        check(random_game(n, it % 2));
    }

    /// A line of 1e6 states, state v moves only to v - 1: the result alternates and moves[v] = v
    int n = 1000000;
    vector<vector<int>> line(n);
    for (int v = 1; v < n; v++) line[v].push_back(v - 1);
    GraphGame chain(line);
    for (int v = 0; v < n; v++){
        assert(chain.result[v] == (v % 2 ? GraphGame::WIN : GraphGame::LOSE));
        assert(chain.moves[v] == v);
    }

    stress::expect_assert("misc/graph_game.cpp", []{ GraphGame bad({{1}, {2}}); });

    return 0;
}
