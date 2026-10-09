/***
 *
 * Knight Distance
 * Fewest knight moves between two squares of an infinite chessboard (every quadrant usable)
 *
 * Complexity: O(1)
 *
 * knight_distance(x, y): moves from (0, 0) to (x, y), any long long coordinates with |x|, |y| <= 1e18
 * For two squares use knight_distance(x2 - x1, y2 - y1)
 *
 * By symmetry only x >= y >= 0 matters. Far from the diagonal (2y < x) the knight walks along x with
 * (2, +-1) steps, near the diagonal it mixes (2, 1) and (1, 2) steps. (1, 0) and (2, 2) are the two exceptions
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long knight_distance(long long x, long long y){
    x = llabs(x), y = llabs(y);
    if (x < y) swap(x, y);
    if (x == 1 && y == 0) return 3;
    if (x == 2 && y == 2) return 4;
    if (y == 0 || 2 * y < x){
        long long c = y & 1, a = x - 2 * c, b = a & 3;
        return (a - b) / 2 + b + c;
    }
    long long d = x - (x - y) / 2, c = (x - y) & 1, z = d % 3 != 0;
    return d / 3 * 2 + c + z * 2 * (1 - c);
}

int main(){
    assert(knight_distance(0, 0) == 0);
    assert(knight_distance(1, 2) == 1 && knight_distance(-2, 1) == 1);
    assert(knight_distance(1, 0) == 3 && knight_distance(0, -1) == 3);
    assert(knight_distance(1, 1) == 2);
    assert(knight_distance(2, 2) == 4);
    assert(knight_distance(3, 3) == 2);
    assert(knight_distance(4, 0) == 2);
    assert(knight_distance(8, 0) == 4);
    assert(knight_distance(7, 7) == 6);
    return 0;
}
