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

    /// Huge coordinates near the diagonal: parity, the lower bound max(x / 2, (x + y) / 3) rounded up, and
    /// a second closed form (delta - 2 floor((delta - y) / k), k = 3 or 4) derived independently of the library's
    auto floor_div = [](long long a, long long b){ return a / b - (a % b != 0 && (a < 0) != (b < 0)); };
    auto reference = [&](long long x, long long y){
        x = llabs(x), y = llabs(y);
        if (x < y) swap(x, y);
        if (x == 1 && y == 0) return 3LL;
        if (x == 2 && y == 2) return 4LL;
        long long delta = x - y;
        return delta - 2 * floor_div(delta - y, y > delta ? 3 : 4);
    };
    for (int x = -CHECK; x <= CHECK; x++){
        for (int y = -CHECK; y <= CHECK; y++) assert(reference(x, y) == dist[(x + R) * W + (y + R)]);
    }
    for (long long it = 0; it < stress::scaled(100000); it++){
        long long x = stress::rand_int(0, 1000000000000000000LL), y = max(0LL, x - stress::rand_int(0, it % 2 ? 1000 : x));
        long long d = knight_distance(x, y), bound = max((x + 1) / 2, (x + y + 2) / 3);
        assert(d == reference(x, y) && d % 2 == (x + y) % 2 && bound <= d && d <= bound + 3);
        assert(d == knight_distance(-x, -y) && d == knight_distance(-y, x));
    }
    return 0;
}
