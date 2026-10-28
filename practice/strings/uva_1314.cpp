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
    /// Two candidate starts i != j, the one that loses a comparison skips past everything it just matched
    /// The smallest least start is never skipped and i starts at 0, so i ends on it
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

    return i;
}

int minimum_rotation(const char* s){
    return minimum_rotation(string(s));
}

int main(){
    static char buf[100005];
    int t, len;

    scanf("%d", &t);
    while (t--){
        scanf("%d %s", &len, buf);
        printf("%d\n", minimum_rotation(string(buf)));
    }

    return 0;
}
