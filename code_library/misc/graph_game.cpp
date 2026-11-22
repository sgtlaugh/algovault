/***
 *
 * Graph Game
 * Win, lose or draw for the player to move in every state of a two-player game whose moves may form cycles
 *
 * Complexity: O(n + m) for n states and m moves
 *
 * GraphGame g(adj): solves the game where the player to move at u picks a state in adj[u], or loses when it is empty
 * g.result[u]: GraphGame::WIN, LOSE or DRAW for the player to move at u
 * g.moves[u]: the number of moves left under optimal play, -1 for a draw
 *
 * Notes:
 *   - Optimal play: the winner ends the game as fast as it can, the loser drags it out as long as it can
 *   - DRAW: neither player can force a win, so the game goes on forever
 *   - Self loops and parallel moves are allowed, every v in adj[u] must lie in [0, n)
 *   - Iterative retrograde analysis from the states without moves, so 1e6 states are safe
 *   - A player who cannot move wins instead: add a state without moves and a move to it from every such state
 *   - With that extra state the results stay right and every finite moves[u] is 1 more than in the real game
 *   - A position where it matters whose turn it is becomes two states, one per player to move
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GraphGame{
    enum Result{LOSE, DRAW, WIN};

    vector<Result> result;
    vector<int> moves;

    GraphGame(const vector<vector<int>>& adj) : result(adj.size(), DRAW), moves(adj.size(), -1){
        int n = adj.size();
        vector<vector<int>> radj(n);
        vector<int> undecided(n), queue;
        for (int u = 0; u < n; u++){
            undecided[u] = adj[u].size();
            for (int v : adj[u]){
                assert(0 <= v && v < n);
                radj[v].push_back(u);
            }
            if (!undecided[u]) result[u] = LOSE, moves[u] = 0, queue.push_back(u);
        }

        /// queue is sorted by moves, so u's first LOSE successor is its closest and its last WIN one its farthest
        for (int i = 0; i < (int)queue.size(); i++){
            int v = queue[i];
            for (int u : radj[v]){
                if (result[u] != DRAW) continue;
                if (result[v] == LOSE || --undecided[u] == 0){
                    result[u] = result[v] == LOSE ? WIN : LOSE;
                    moves[u] = moves[v] + 1;
                    queue.push_back(u);
                }
            }
        }
    }
};

int main(){
    /***
     * 0 -> 1 -> 2, 2 has no move; 3 -> 3 forever; 4 -> 0 or 3; 5 -> 1 or 2; 6 -> 1 or 4; 7 -> 8, 8 -> 7 or 2
     * 2 loses at once, 1 wins in 1, 0 loses in 2, 4 wins by moving to 0 in 3, 5 wins in 1 by moving to 2
     * 6 loses either way and drags it out through 4: 1 + 3 = 4; 8 escapes the cycle to win in 1, so 7 loses in 2
    ***/
    vector<vector<int>> adj = {{1}, {2}, {}, {3}, {0, 3}, {1, 2}, {1, 4}, {8}, {7, 2}};
    GraphGame g(adj);
    assert(g.result[2] == GraphGame::LOSE);
    assert(g.moves[2] == 0);
    assert(g.result[0] == GraphGame::LOSE);
    assert(g.moves[0] == 2);
    assert(g.result[3] == GraphGame::DRAW);
    assert(g.moves[3] == -1);
    assert(g.result[4] == GraphGame::WIN);
    assert(g.moves[4] == 3);
    assert(g.moves[5] == 1);
    assert(g.result[6] == GraphGame::LOSE);
    assert(g.moves[6] == 4);
    assert(g.result[7] == GraphGame::LOSE);
    assert(g.moves[7] == 2);
    assert(g.result[8] == GraphGame::WIN);

    /// Two states chasing each other forever: both draw
    GraphGame ring({{1}, {0}});
    assert(ring.result[0] == GraphGame::DRAW);
    assert(ring.result[1] == GraphGame::DRAW);

    return 0;
}
