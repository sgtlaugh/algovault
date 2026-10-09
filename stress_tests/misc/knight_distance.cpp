#include "../common.h"

#define main library_main
#include "../../code_library/misc/knight_distance.cpp"
#undef main

/// BFS from the origin over a board wide enough that no shortest path to the checked squares leaves it
int main(){
    const int R = 140, CHECK = 90, W = 2 * R + 1;
    vector<int> dist(W * W, -1);
    vector<int> queue = {R * W + R};
    dist[R * W + R] = 0;
    const int dx[] = {1, 2, 2, 1, -1, -2, -2, -1}, dy[] = {2, 1, -1, -2, -2, -1, 1, 2};
    for (size_t i = 0; i < queue.size(); i++){
        int x = queue[i] / W, y = queue[i] % W;
        for (int k = 0; k < 8; k++){
            int nx = x + dx[k], ny = y + dy[k];
            if (nx < 0 || ny < 0 || nx >= W || ny >= W || dist[nx * W + ny] != -1) continue;
            dist[nx * W + ny] = dist[queue[i]] + 1;
            queue.push_back(nx * W + ny);
        }
    }
    for (int x = -CHECK; x <= CHECK; x++){
        for (int y = -CHECK; y <= CHECK; y++) assert(knight_distance(x, y) == dist[(x + R) * W + (y + R)]);
    }

    /// Huge coordinates: the formula must agree with its own translation by whole (4, 2) blocks of two moves
    for (long long it = 0; it < stress::scaled(100000); it++){
        long long x = stress::rand_int(0, 1000000000000000000LL), y = stress::rand_int(-x / 3 - 1, x / 3 + 1);
        long long d = knight_distance(x, y);
        assert(d >= 0 && d == knight_distance(-x, y) && d == knight_distance(y, x) && d == knight_distance(x, -y));
        if (x > 100 && 2 * llabs(y) + 10 < x) assert(knight_distance(x + 4, y) == d + 2);
    }
    return 0;
}
