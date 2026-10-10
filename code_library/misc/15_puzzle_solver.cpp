/***
 *
 * 15 Puzzle Solver
 * IDA* with the Manhattan distance heuristic, finds an optimal move sequence for a 4x4 board
 *
 * Complexity: O(3^d) for an optimal solution of d moves, the threshold strictly increases up to d
 *   and an iteration with threshold t visits O(3^t) states since the blank never steps straight back, O(d) memory
 *
 * FifteenPuzzle p(tiles): tiles holds the 16 cells row by row, 0 is the blank
 * p.is_solvable(): true if the board can reach the solved state 1 2 ... 15 0
 * p.solve(): {move count, moves} of an optimal solution, {-1, ""} if the board is not solvable
 *   each move is the direction the blank slides: L, R, U or D
 *
 * Example:
 *   FifteenPuzzle p({2, 3, 4, 0, 1, 5, 7, 8, 9, 6, 10, 12, 13, 14, 11, 15});
 *   auto [moves, path] = p.solve();  // moves = 9, path = "LLLDRDRDR"
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FifteenPuzzle{
    static constexpr char dir[] = "LRUD";
    static constexpr int dx[] = {0, 0, -1, 1};
    static constexpr int dy[] = {-1, 1, 0, 0};

    int ar[4][4], blank_x = 0, blank_y = 0;  /// blank stored as 16 so (ar - 1) is the tile's home cell
    string path;

    FifteenPuzzle(const vector<int>& tiles){
        assert(tiles.size() == 16);
        for (int i = 0; i < 16; i++){
            ar[i >> 2][i & 3] = tiles[i] ? tiles[i] : 16;
            if (!tiles[i]) blank_x = i >> 2, blank_y = i & 3;
        }
    }

    bool is_solvable() const{
        int i, j, r = 0, counter = 0;

        for (i = 0; i < 16; i++){
            if (ar[i >> 2][i & 3] == 16) r = (i >> 2);
            else{
                for (j = i + 1; j < 16; j++){
                    if (ar[j >> 2][j & 3] < ar[i >> 2][i & 3]) counter++;
                }
            }
        }

        return ((counter + r) & 1);
    }

    pair<int, string> solve(){
        if (!is_solvable()) return {-1, ""};

        int lim = heuristic();
        while (true){
            path.assign(lim, 0);
            int res = ida_star(blank_x, blank_y, blank_x, blank_y, 0, lim, heuristic());
            if (res <= lim) return {res, path};
            lim = res;
        }
    }

    int get_cost(int i, int j) const{
        return abs(((ar[i][j] - 1) >> 2) - i) + abs(((ar[i][j] - 1) & 3) - j);
    }

    int heuristic() const{
        int i, j, manhattan_distance = 0;

        for (i = 0; i < 4; i++){
            for (j = 0; j < 4; j++){
                if (ar[i][j] & 15) manhattan_distance += get_cost(i, j);
            }
        }

        return manhattan_distance;
    }

    /// Returns g once solved, which is <= lim and makes every caller return at once, else the next threshold
    int ida_star(int bx, int by, int lx, int ly, int g, int lim, int h){
        if (!h){
            path.resize(g);
            return g;
        }

        int f = g + h;
        if (f > lim) return f;

        int i, k, l, nh, r, res = 1 << 20;
        for (i = 0; i < 4; i++){
            k = bx + dx[i], l = by + dy[i];
            if (k >= 0 && k < 4 && l >= 0 && l < 4 && !(k == lx && l == ly)){
                nh = h;
                nh -= get_cost(k, l);
                swap(ar[bx][by], ar[k][l]);
                nh += get_cost(bx, by);

                path[g] = dir[i];
                r = ida_star(k, l, bx, by, g + 1, lim, nh);
                swap(ar[bx][by], ar[k][l]);
                if (r < res) res = r;
                if (r <= lim) return r;
            }
        }

        return res;
    }
};

int main(){
    assert((FifteenPuzzle({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0}).solve() == pair<int, string>{0, ""}));
    assert((FifteenPuzzle({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 15}).solve() == pair<int, string>{1, "R"}));
    assert((FifteenPuzzle({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 13, 14, 15, 12}).solve() == pair<int, string>{1, "D"}));
    assert((FifteenPuzzle({2, 3, 4, 0, 1, 5, 7, 8, 9, 6, 10, 12, 13, 14, 11, 15}).solve() == pair<int, string>{9, "LLLDRDRDR"}));

    FifteenPuzzle swapped({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 14, 0});
    assert(!swapped.is_solvable());
    assert((swapped.solve() == pair<int, string>{-1, ""}));

    FifteenPuzzle unsolvable({13, 1, 2, 4, 5, 0, 3, 7, 9, 6, 10, 12, 15, 8, 11, 14});
    assert(!unsolvable.is_solvable());
    assert((unsolvable.solve() == pair<int, string>{-1, ""}));

    FifteenPuzzle hard({6, 2, 8, 4, 12, 14, 1, 10, 13, 15, 3, 9, 11, 0, 5, 7});
    assert(hard.is_solvable());
    assert((hard.solve() == pair<int, string>{48, "RULDRURDLLLURURULDRDDLUURDLULURRDRDLULLDRRURDLDR"}));
    assert((hard.solve().first == 48));

    return 0;
}
