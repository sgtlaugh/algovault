/***
 *
 * Z Algorithm
 * Computes Z array where Z[i] = length of longest common prefix of s[i..] and s
 *
 * Complexity: O(n)
 *
 * Algorithm maintains window [l, r] where s[l..r] matches prefix s[0..r-l]
 * For each position i, reuses previously computed values when inside window
 * Z[i] = length of longest prefix of s starting at i that matches s[0..]
 *
 * Applications:
 *   - Pattern matching: concatenate "pattern$text", Z[i] == |pattern| means match at (i - |pattern| - 1)
 *   - Find period of string
 *   - String borders/prefixes analysis
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Container>
vector<int> z_function(const Container& c){
    int n = c.size();
    vector<int> z(n);

    int l = 0, r = 0;
    for (int i = 1; i < n; i++){
        int val = z[i - l];
        if (i + val >= r){
            val = max(r - i, 0);
            while (i + val < n && c[val] == c[i + val]){
                val++;
            }
            l = i;
            r = i + val;
        }
        z[i] = val;
    }

    z[0] = n;
    return z;
}

vector<int> z_function(const char* s){
    return z_function(string(s));
}

int main(){
    assert((z_function("a") == vector<int>{1}));
    assert((z_function("aaaa") == vector<int>{4, 3, 2, 1}));
    assert((z_function("aaabaab") == vector<int>{7, 2, 1, 0, 2, 1, 0}));
    assert((z_function("abababab") == vector<int>{8, 0, 6, 0, 4, 0, 2, 0}));
    assert((z_function("abracadabra") == vector<int>{11, 0, 0, 1, 0, 1, 0, 4, 0, 0, 1}));

    assert((z_function(vector<int>{1, 2, 1, 2, 3}) == vector<int>{5, 0, 2, 0, 0}));

    return 0;
}
