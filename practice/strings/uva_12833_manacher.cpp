/***
 *
 * Manacher's Algorithm
 * Length of the longest palindrome around every center of a string or vector
 *
 * Complexity: O(n) time and memory
 *
 * pal = manacher(s) has 2n - 1 entries, empty for an empty s, 0-based positions, ranges [l, r] inclusive
 *   i even: pal[i] = length of the longest odd palindrome centered at s[i / 2]
 *   i odd:  pal[i] = length of the longest even palindrome centered between s[i / 2] and s[i / 2 + 1], 0 if none
 * s[l..r] is a palindrome iff pal[l + r] >= r - l + 1
 * The palindrome at center i starts at (i - pal[i] + 1) / 2, the longest one has length *max_element(pal)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Container>
vector<int> manacher(const Container& str){
    if (str.empty()) return {};
    int i, j, k, l = str.size(), n = l << 1;
    vector<int> pal(n);

    for (i = 0, j = 0, k = 0; i < n; j = max(0, j - k), i += k){
        while (j <= i && (i + j + 1) < n && str[(i - j) >> 1] == str[(i + j + 1) >> 1]) j++;
        for (k = 1, pal[i] = j; k <= i && k <= pal[i] && (pal[i] - k) != pal[i - k]; k++){
            pal[i + k] = min(pal[i - k], pal[i] - k);
        }
    }

    pal.pop_back();
    return pal;
}

vector<int> manacher(const char* s){
    return manacher(string(s));
}

int main(){
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, q;
        static char buf[100010];
        scanf("%d %s", &n, buf);
        string s = buf;

        vector<array<int, 26>> pre(n + 1);
        pre[0].fill(0);
        for (int i = 0; i < n; i++) pre[i + 1] = pre[i], pre[i + 1][s[i] - 'a']++;

        vector<vector<int>> odd(26, vector<int>(n + 2, 0)), even(26, vector<int>(n + 2, 0));
        vector<int> pal = manacher(s);
        for (int i = 0; i < 2 * n - 1; i++){
            if (!pal[i]) continue;
            int l = (i - pal[i] + 1) / 2, mid = i / 2;
            for (int c = 0; c < 26; c++){
                int x = pre[mid + 1][c] - pre[l][c];
                if (i % 2 == 0 && s[mid] - 'a' == c) odd[c][x]++;
                else even[c][x]++;
            }
        }
        for (int c = 0; c < 26; c++){
            for (int x = n; x >= 1; x--) odd[c][x] += odd[c][x + 1], even[c][x] += even[c][x + 1];
        }

        static char qs[100010];
        scanf("%d %s", &q, qs);
        printf("Case %d:\n", tc);
        for (int i = 0; i < q; i++){
            int x, c = qs[i] - 'a';
            scanf("%d", &x);
            int res = 0;
            if (x > 0 && x <= n) res = x % 2 ? odd[c][(x + 1) / 2] : even[c][x / 2];
            printf("%d\n", res);
        }
    }

    return 0;
}
