/***
 *
 * Z Algorithm
 *
 * z_function(s): z[i] = length of the longest common prefix of s[i..] and s, z[0] = n, empty for n = 0
 *     s is a string, C string, vector or any container with size() and [] whose elements compare with ==
 *
 * Complexity: O(n)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Container>
vector<int> z_function(const Container& c){
    int n = c.size();
    vector<int> z(n);
    if (n == 0) return z;

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
    int t;
    cin >> t;
    for (int tc = 0; tc < t; tc++){
        string s;
        cin >> s;
        int n = s.size(), p = 1;
        vector<int> z = z_function(s);
        /// s repeats its first p characters iff p divides n and s[p..] matches a prefix all the way to the end
        while (p < n && (n % p || z[p] != n - p)) p++;
        if (tc) puts("");
        printf("%d\n", p);
    }
    return 0;
}
