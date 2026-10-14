/***
 *
 * Minimum Rotation
 * Start of the lexicographically smallest rotation of a sequence
 *
 * Complexity: O(n), O(1) extra memory
 *
 * minimum_rotation(s): smallest index i such that s[i..] + s[..i) is the least rotation, 0 for an empty s
 * Works on strings, vectors, C strings (literals, char arrays, char pointers) and any indexable container whose elements compare with < and ==
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

/// Without the char* overload a non-const char[] or char* skips the const conversion by binding to the template, which has no size()
int minimum_rotation(const char* s){
    return minimum_rotation(string_view(s));
}

int minimum_rotation(char* s){
    return minimum_rotation(string_view(s));
}

int main(){
    /// The rotations of "bca" are bca, cab, abc, so the least starts at 2
    assert(minimum_rotation(string("bca")) == 2);

    /// "abab" starts at 1 and at 3, the smaller index wins
    assert(minimum_rotation(string("baba")) == 1);

    /// C strings: a literal, then a char array longer than its contents
    assert(minimum_rotation("dcba") == 3);

    char word[16] = "cab";
    assert(minimum_rotation(word) == 1);

    /// Any indexable container, the least rotation is 1 1 3 1 2
    assert(minimum_rotation(vector<int>{3, 1, 2, 1, 1}) == 3);

    return 0;
}
