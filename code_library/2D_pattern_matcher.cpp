#include <bits/stdtr1c++.h>

#define MAX 1010


using namespace std;

/***
 * 2D Pattern Matcher, n x m text, r x c pattern
 * Processes and stores the output in the occur[][]
 * occur[i][j] = true for every (k, l), where k >= 0 && k < r && l >= 0 && l < c, text[i + k][j + l] == pattern[k][l]
 *
 * For example:
 *
 * text:
 * abab
 * baba
 * abaa
 *
 * pattern:
 * ba
 * ab
 *
 * occur matrix:
 * 0100
 * 1000
 * 0000

***/

namespace pm{
    const unsigned long long base = 1995433697ULL;
    const unsigned long long base2 = 1000000007ULL;

    int n, m, r, c, match;
    bool occur[MAX][MAX];

    char str[MAX][MAX], pattern[MAX][MAX];
    unsigned long long pattern_hash, P[MAX], F[MAX], ar[MAX];

	/// Call init once with the text and pattern
    void init(int a, int b, char text[MAX][MAX], int u, int v, char pat[MAX][MAX]){
        n = a, m = b, r = u, c = v;
        for (int i = 0; i < n; i++) text[i][m] = 0, strcpy(str[i], text[i]);
        for (int i = 0; i < r; i++) pat[i][c] = 0, strcpy(pattern[i], pat[i]);

        P[0] = F[0] = 1;
        for (int i = 1; i < MAX; i++){
            P[i] = P[i - 1] * base;
            F[i] = F[i - 1] * base2;
        }

        pattern_hash = 0;
        for (int i = 0; i < r; i++){
            unsigned long long res = 0;
            for (int j = 0; j < c; j++) res = (res * base2) + pattern[i][j];
            pattern_hash = (pattern_hash * base) + res;
        }

        for (int i = 0; i < n; i++){
            unsigned long long res = 0;
            for (int j = 0; j < c; j++) res = (res * base2) + str[i][j];
            ar[i] = res;
        }
    }

    void roll(int col){
        unsigned long long x = 0, y;
        for (int i = 0; i < r; i++) x = (x * base) + ar[i];

        for (int i = 0; i < n; i++){
            if ((i + r) > n) break;
            if (x == pattern_hash) match++, occur[i][col] = true;

            y = (x - (P[r - 1] * ar[i]));
            x = (y * base) + ar[i + r];
        }
    }

	/// Call solve once to populate the occur matrix
	/// Can configure to return the matched positions in a vector of pais if required
    void solve(){
        match = 0;
        memset(occur, 0, sizeof(occur));

        for (int j = 0; j < m; j++){
            if ((j + c) > m) break;

            roll(j);
            for (int i = 0; i < n; i++){
                unsigned long long x = ar[i];
                x = x - (F[c - 1] * str[i][j]);
                ar[i] = (x * base2) + str[i][j + c];
            }
        }
    }
}

int n, m, r, c;
char text[MAX][MAX], pattern[MAX][MAX];

int main(){
	/***

	3 4
	abab
	baba
	abaa

	2 2
	ba
	ab

	***/


	int i, j;
	scanf("%d %d", &n, &m);
	for (i = 0; i < n; i++) scanf("%s", text[i]);

	scanf("%d %d", &r, &c);
	for (i = 0; i < r; i++) scanf("%s", pattern[i]);

	pm::init(n, m, text, r, c, pattern);
	pm::solve();

	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			printf("%d", pm::occur[i][j]);
		}
		puts("");
	}

	return 0;
}
