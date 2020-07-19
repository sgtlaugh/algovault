#include <bits/stdc++.h>

using namespace std;

/***
 *
 * Manacher's algorithm to generate longest palindromic substrings for all centers
 * For a string of length n, returns a vector of size 2 * n - 1
 * When i is even, pal[i] = largest palindromic substring centered from str[i / 2]
 * When i is odd,  pal[i] = largest palindromic substring centered between str[i / 2] and str[i / 2] + 1
 *
***/

vector <int> manacher(const string& str){
    int i, j, k, l = str.size(), n = l << 1;
    vector <int> pal(n);

    for (i = 0, j = 0, k = 0; i < n; j = max(0, j - k), i += k){
        while (j <= i && (i + j + 1) < n && str[(i - j) >> 1] == str[(i + j + 1) >> 1]) j++;
        for (k = 1, pal[i] = j; k <= i && k <= pal[i] && (pal[i] - k) != pal[i - k]; k++){
            pal[i + k] = min(pal[i - k], pal[i] - k);
        }
    }

    pal.pop_back();
    return pal;
}

int main(){
    auto v = manacher("abababbaa");
    for (auto x: v) printf("%d ", x);  /// 1 0 3 0 5 0 5 0 3 0 1 4 1 0 1 2 1
    puts("");

    return 0;
}
