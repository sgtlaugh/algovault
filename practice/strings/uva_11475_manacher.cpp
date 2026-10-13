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
    static char buf[100010];
    while (scanf("%s", buf) == 1){
        string s = buf;
        int n = s.size(), i = 0;
        vector<int> pal = manacher(s);

        /// Mirror the part before the longest palindromic suffix
        while (pal[i + n - 1] < n - i) i++;
        puts((s + string(s.rend() - i, s.rend())).c_str());
    }
    return 0;
}
