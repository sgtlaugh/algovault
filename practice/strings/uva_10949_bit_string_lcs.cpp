/***
 *
 * Bit-String LCS
 * Calculates the longest common subsequence of two strings with bitset
 * L. Allison and T. I. Dix, "A bit-string longest-common-subsequence algorithm", Information Processing Letters 23, pp. 305-310, 1986, doi:10.1016/0020-0190(86)90091-8
 *
 * lcs(A, B): the LCS length of the null-terminated strings A and B (const char*), any byte values except 0
 *
 * Complexity: O(N * M / 64)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Hyyro's row update V = (V + U) | (V - U) with U = V & match, zeros of V count the LCS
/// U is a submask of V so V - U never borrows, only the addition carries into the next 64-bit block
int lcs(const char* A, const char* B){
    int i, j, k, n = strlen(A), m = strlen(B), res = 0;
    unsigned long long u, v, sum, mask[256];
    vector<char> carry(n, 0);

    for (i = 0; i * 64 < m; i++){
        memset(mask, 0, sizeof(mask));
        for (k = 0; k < 64 && i * 64 + k < m; k++){
            mask[(unsigned char)B[i * 64 + k]] |= (1ULL << k);
        }

        /// Bits past m stay set, since mask is 0 there and V - U keeps them
        for (j = 0, v = ~0ULL; j < n; j++){
            u = v & mask[(unsigned char)A[j]];
            bool c = __builtin_uaddll_overflow(v, u, &sum);
            c |= __builtin_uaddll_overflow(sum, carry[j], &sum);
            carry[j] = c;
            v = sum | (v - u);
        }

        res += __builtin_popcountll(~v);
    }

    return res;
}

char grid[25][25];

string walk(){
    static char line[20010];
    int n, x, y;

    scanf("%d %d %d", &n, &x, &y);
    fgets(line, sizeof(line), stdin);
    if (!fgets(line, sizeof(line), stdin)) line[0] = 0;

    string s(1, grid[--x][--y]);
    for (int i = 0; i < n; i++){
        if (line[i] == 'N') x--;
        if (line[i] == 'S') x++;
        if (line[i] == 'E') y++;
        if (line[i] == 'W') y--;
        s += grid[x][y];
    }

    return s;
}

int main(){
    int t, h, w;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        scanf("%d %d", &h, &w);
        for (int i = 0; i < h; i++) scanf("%s", grid[i]);

        string a = walk(), b = walk();
        int l = lcs(a.c_str(), b.c_str());
        printf("Case %d: %d %d\n", tc, (int)a.size() - l, (int)b.size() - l);
    }

    return 0;
}
