/***
 *
 * Manacher's Algorithm
 * Length of the longest palindrome around every center of a string or vector
 *
 * Complexity: O(n)
 *
 * pal = manacher(s): 2n - 1 entries, empty for an empty s
 * pal[i] for even i: length of the longest odd palindrome centered at s[i / 2]
 * pal[i] for odd i: length of the longest even palindrome centered between s[i / 2] and s[i / 2 + 1], 0 if none
 *
 * Notes:
 *   - s[l..r], l and r inclusive, is a palindrome iff pal[l + r] >= r - l + 1
 *   - The palindrome at center i starts at (i - pal[i] + 1) / 2
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
    assert(manacher("abababbaa") == vector<int>({1, 0, 3, 0, 5, 0, 5, 0, 3, 0, 1, 4, 1, 0, 1, 2, 1}));
    assert(manacher("") == vector<int>{});
    assert(manacher("z") == vector<int>{1});
    assert(manacher(vector<int>{4, 9, 9, 4}) == vector<int>({1, 0, 1, 4, 1, 0, 1}));

    return 0;
}
