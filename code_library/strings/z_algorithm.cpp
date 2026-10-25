/***
 *
 * Z Algorithm
 * Longest common prefix of every suffix with the whole sequence
 *
 * Complexity: O(n)
 *
 * z_function(s): z[i] = length of the longest common prefix of s[i..] and s, z[0] = n, empty for n = 0
 *
 * Notes:
 *   - s is a string, C string, vector or any container with size() and [] whose elements compare with ==
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
    assert((z_function("a") == vector<int>{1}));
    assert((z_function("aaaa") == vector<int>{4, 3, 2, 1}));
    assert((z_function("aaabaab") == vector<int>{7, 2, 1, 0, 2, 1, 0}));
    assert((z_function("abababab") == vector<int>{8, 0, 6, 0, 4, 0, 2, 0}));
    assert((z_function("abracadabra") == vector<int>{11, 0, 0, 1, 0, 1, 0, 4, 0, 0, 1}));

    assert((z_function(vector<int>{1, 2, 1, 2, 3}) == vector<int>{5, 0, 2, 0, 0}));

    return 0;
}
