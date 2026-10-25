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
    static char buf[1000010];
    for (int tc = 1; tc <= t; tc++){
        int x, y;
        scanf("%d %d %s", &x, &y, buf);
        string s = buf;
        int n = s.size();
        vector<int> pal = manacher(s);

        vector<int> first(n + 1, -1);
        for (int st = 1, covered = 0; st < n; st++){
            int reach = st + pal[2 * st - 1] / 2;
            for (int e = max(st, covered) + 1; e <= reach; e++) first[e] = st;
            covered = max(covered, reach);
        }

        vector<long long> best(n + 1, 0);
        for (int e = 1; e <= n; e++){
            best[e] = best[e - 1] + x;
            if (first[e] != -1) best[e] = min(best[e], best[first[e]] + y);
        }
        printf("Case %d: %lld\n", tc, best[n]);
    }
    return 0;
}
