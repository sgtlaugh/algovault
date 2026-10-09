/***
 *
 * Minimum Rotation
 * Start of the lexicographically smallest rotation of a sequence
 *
 * Complexity: O(n), O(1) extra memory
 *
 * minimum_rotation(s): smallest index i such that s[i..] + s[..i) is the least rotation, 0 for an empty s
 * Works on strings, vectors and any indexable container whose elements compare with <
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Container>
int minimum_rotation(const Container& s){
    int n = s.size(), i = 0, j = 1, k = 0;
    /// Two candidate starts i < j, the one that loses a comparison skips past everything it just matched
    while (i < n && j < n && k < n){
        auto a = s[(i + k) % n], b = s[(j + k) % n];
        if (a == b){
            k++;
            continue;
        }
        if (b < a) i += k + 1;
        else j += k + 1;
        if (i == j) j++;
        k = 0;
    }
    return n == 0 ? 0 : min(i, j);
}

int minimum_rotation(const char* s){
    return minimum_rotation(string(s));
}

int main(){
    assert(minimum_rotation(string("bca")) == 2);
    assert(minimum_rotation(string("abc")) == 0);
    assert(minimum_rotation(string("cab")) == 1);
    assert(minimum_rotation(string("abab")) == 0);
    assert(minimum_rotation(string("baba")) == 1);
    assert(minimum_rotation(string("aaaa")) == 0);
    assert(minimum_rotation(string("")) == 0);
    assert(minimum_rotation(string("z")) == 0);
    assert(minimum_rotation("dcba") == 3);
    assert(minimum_rotation(vector<int>{3, 1, 2, 1, 1}) == 3);
    return 0;
}
