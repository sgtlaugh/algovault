#include "../common.h"

#define main library_main
#include "../../code_library/misc/15_puzzle_solver.cpp"
#undef main

const string IN = "/tmp/algovault_15_puzzle_in.txt", OUT = "/tmp/algovault_15_puzzle_out.txt";
const int bx_step[] = {0, 0, -1, 1}, by_step[] = {-1, 1, 0, 0};  /// L R U D move the blank, as in the solver's output
typedef array<int, 16> Board;  /// row major, 0 is the blank

int blank(const Board& b){
    return find(b.begin(), b.end(), 0) - b.begin();
}

bool slide(Board& b, int d){
    int p = blank(b), r = p / 4 + bx_step[d], c = p % 4 + by_step[d];
    if (r < 0 || r > 3 || c < 0 || c > 3) return false;
    swap(b[p], b[r * 4 + c]);
    return true;
}

Board goal(){
    Board b;
    for (int i = 0; i < 16; i++) b[i] = (i + 1) % 16;
    return b;
}

int bfs_distance(const Board& start){
    map<Board, int> dist = {{start, 0}};
    queue<Board> q;
    q.push(start);

    while (true){
        Board b = q.front();
        q.pop();
        if (b == goal()) return dist[b];
        for (int d = 0; d < 4; d++){
            Board n = b;
            if (slide(n, d) && !dist.count(n)) dist[n] = dist[b] + 1, q.push(n);
        }
    }
}

/// Solvable exactly when the permutation parity (blank counted as 16) matches the parity of the blank's distance from its home
bool solvable_oracle(const Board& b){
    int perm[16], parity = 0;
    for (int i = 0; i < 16; i++) perm[i] = b[i] ? b[i] - 1 : 15;
    for (int i = 0; i < 16; i++) for (int j = i + 1; j < 16; j++) parity ^= perm[i] > perm[j];

    int p = blank(b);
    return parity == ((3 - p / 4 + 3 - p % 4) & 1);
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        Board b;
        iota(b.begin(), b.end(), 0);
        shuffle(b.begin(), b.end(), stress::rng());
        for (int i = 0; i < 16; i++) ar[i / 4][i % 4] = b[i] ? b[i] : 16;
        assert(is_solvable() == solvable_oracle(b));
    }

    for (long long it = 0; it < stress::scaled(10); it++){
        vector<Board> boards;
        vector<int> walks;
        vector<bool> solvable;
        for (int k = 0; k < 20; k++){
            Board b = goal();
            int walk = stress::rand_int(0, k % 2 ? 16 : 30);
            for (int s = 0; s < walk; ){
                if (slide(b, stress::rand_int(0, 3))) s++;
            }
            bool ok = k % 5 != 4;
            if (!ok){
                int x = stress::rand_int(1, 15), y = x % 15 + 1;  /// a transposition of two tiles flips solvability
                swap(*find(b.begin(), b.end(), x), *find(b.begin(), b.end(), y));
            }
            boards.push_back(b), walks.push_back(walk), solvable.push_back(ok);
        }

        string input = to_string(boards.size()) + "\n";
        for (auto& b : boards) for (int i = 0; i < 16; i++) input += to_string(b[i]) + (i % 4 == 3 ? "\n" : " ");
        ofstream(IN) << input;
        assert(freopen(IN.c_str(), "r", stdin) && freopen(OUT.c_str(), "w", stdout));
        library_main();
        fflush(stdout);

        ifstream out(OUT);
        string line;
        for (size_t k = 0; k < boards.size(); k++){
            while (getline(out, line) && line.empty()) {}
            if (!solvable[k]){
                assert(line == "This puzzle is not solvable.");
                continue;
            }
            int moves = -1;
            assert(sscanf(line.c_str(), "Puzzle can be solved in %d moves", &moves) == 1);
            string path;
            if (moves) getline(out, path);
            assert((int)path.size() == moves && moves <= walks[k]);

            Board b = boards[k];
            for (char ch : path) assert(slide(b, strchr("LRUD", ch) - "LRUD"));
            assert(b == goal());
            if (walks[k] <= 16) assert(moves == bfs_distance(boards[k]));  /// IDA* overshooting its bound returns longer paths
        }
    }

    return 0;
}
