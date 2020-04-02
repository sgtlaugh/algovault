#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;


const char dir[] = "LRUD";
const int dx[] = {0, 0, -1, 1};
const int dy[] = {-1, 1, 0, 0};

char str[60];
int found, len[4] = {0}, idx[4][4], ar[4][4];


int get_cost(int i, int j){
	return abs(((ar[i][j] - 1) >> 2) - i) + abs( ((ar[i][j] - 1) & 3) - j);
}

int heuristic(int bx, int by){
    int i, j, res, manhattan_distance = 0;

    for (i = 0; i < 4; i++){
        for (j = 0; j < 4; j++){
            if (ar[i][j] & 15){
                manhattan_distance += get_cost(i, j);
            }
        }
    }

    res = manhattan_distance;
    return res;
}

int ida_star(int bx, int by, int lx, int ly, int g, int lim, int d, int h){
    if (found) return 1 << 30;

    if (!h){
		found = 1;
		str[d] = 0;
		printf("Puzzle can be solved in %d moves\n", g);
		printf("%s\n\n", str);
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

            str[d] = dir[i];
            r = ida_star(k, l, bx, by, g + 1, lim, d + 1, nh);
            swap(ar[bx][by], ar[k][l]);
            if (r < res) res = r;
            if (r <= lim) return r;
        }
    }

    return res;
}

void solve(int bx, int by){
    int res, lim = heuristic(bx, by);

    found = 0;
    while (!found) {
        res = ida_star(bx, by, bx, by, 0, lim, 0, heuristic(bx, by));
        if (res <= lim) return;
        else lim = res;
    }
}

bool is_solvable(int bx, int by){
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

int main(){
    int t, i, j, bx, by;

    scanf("%d", &t);
    while (t--){
        for (i = 0; i < 4; i++){
            for (j = 0; j < 4; j++){
                scanf("%d", &ar[i][j]);
                if (!ar[i][j]){
                    bx = i, by = j;
                    ar[i][j] = 16;
                }
            }
        }

        if (!is_solvable(bx, by)) puts("This puzzle is not solvable.\n");
        else solve(bx, by);
    }
    return 0;
}

/***

Input:


6

2 3 4 0
1 5 7 8
9 6 10 12
13 14 11 15


6 2 8 4
12 14 1 10
13 15 3 9
11 0 5 7

6 8 4 2
12 14 1 10
13 15 3 9
11 0 5 7

13 1 2 4
5 0 3 7
9 6 10 12
15 8 11 14

1 8 14 12
11 2 9 10
0 7 6 5
15 13 4 3

6 3 1 9
15 13 12 11
5 8 10 14
4 0 7 2


Output:


Puzzle can be solved in 9 moves
LLLDRDRDR

Puzzle can be solved in 48 moves
RULDRURDLLLURURULDRDDLUURDLULURRDRDLULLDRRURDLDR

Puzzle can be solved in 50 moves
RRULLDLURURRDLDRULDLURUULDLURRRDLLURRDLLLDRRRULDDR

This puzzle is not solvable.

Puzzle can be solved in 54 moves
RRUULDLDRRDLLURURDDLURRDLLUULDDRURRUULDDRUULDLDLURDRDR

Puzzle can be solved in 60 moves
LURRDLURDRULLUURRDDLUURDDLDLULUURDRDLLDRRRULLURULDLDRDRUURDD



***/
