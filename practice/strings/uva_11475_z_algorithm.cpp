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
    static char buf[100010];

    while (scanf("%s", buf) == 1){
        string s = buf, r(s.rbegin(), s.rend());
        int n = s.size(), i = 0;

        vector<int> z = z_function(r + "#" + s);
        while (z[n + 1 + i] < n - i) i++;
        puts((s + string(s.rend() - i, s.rend())).c_str());
    }

    return 0;
}
