/***
 *
 * Lyndon Factorization (Duval)
 * Splits a sequence into Lyndon words w1 >= w2 >= ... >= wk whose concatenation is the sequence
 *
 * Complexity: O(n), O(1) extra memory besides the output, the const char* overloads copy the input
 *
 * A Lyndon word is strictly smaller than each of its proper suffixes, the factorization is unique
 * lyndon_factorization(s): start index of every factor in increasing order, factor t is s[starts[t], starts[t + 1])
 *     with n closing the last one, empty for an empty s
 * lyndon_min_rotation(s): smallest index i such that s[i..] + s[..i) is the least rotation, 0 for an empty s
 *     Duval on s + s read cyclically, the least rotation starts at the first copy of the last factor group that begins below n
 * Works on strings, vectors and any indexable container whose elements compare with <
 *
 * lyndon_factorization("banana") = {0, 1, 3, 5}       b | an | an | a
 * lyndon_min_rotation("banana") = 5                   abanan
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Container>
vector<int> lyndon_factorization(const Container& s){
    int n = s.size();
    vector<int> starts;

    /// s[i, j) is a power of the Lyndon word s[i, i + j - k) plus a prefix of it, emit the full copies
    for (int i = 0; i < n;){
        int j = i + 1, k = i;
        while (j < n && !(s[j] < s[k])){
            if (s[k] < s[j]) k = i;
            else k++;
            j++;
        }
        for (; i <= k; i += j - k) starts.push_back(i);
    }

    return starts;
}

vector<int> lyndon_factorization(const char* s){
    return lyndon_factorization(string(s));
}

template <typename Container>
int lyndon_min_rotation(const Container& s){
    int n = s.size(), best = 0;

    /// best is the first copy of a repeated factor, so equal least rotations resolve to the smallest index
    for (int i = 0; i < n;){
        best = i;
        int j = i + 1, k = i;
        while (j < 2 * n && !(s[j % n] < s[k % n])){
            if (s[k % n] < s[j % n]) k = i;
            else k++;
            j++;
        }
        while (i <= k) i += j - k;
    }

    return best;
}

int lyndon_min_rotation(const char* s){
    return lyndon_min_rotation(string(s));
}

int main(){
    assert((lyndon_factorization("banana") == vector<int>{0, 1, 3, 5}));
    assert((lyndon_factorization("bara") == vector<int>{0, 1, 3}));
    assert((lyndon_factorization("abracadabra") == vector<int>{0, 7, 10}));
    assert((lyndon_factorization("abab") == vector<int>{0, 2}));
    assert((lyndon_factorization("aab") == vector<int>{0}));
    assert((lyndon_factorization("aaa") == vector<int>{0, 1, 2}));
    assert((lyndon_factorization("dcba") == vector<int>{0, 1, 2, 3}));
    assert((lyndon_factorization("a") == vector<int>{0}));
    assert((lyndon_factorization("") == vector<int>{}));
    assert((lyndon_factorization(vector<int>{3, 1, 2, 1, 1}) == vector<int>{0, 1, 3, 4}));
    assert((lyndon_factorization(vector<int>{-5, 0, -5, 0, -5}) == vector<int>{0, 2, 4}));

    assert(lyndon_min_rotation("banana") == 5);
    assert(lyndon_min_rotation("bca") == 2);
    assert(lyndon_min_rotation("abc") == 0);
    assert(lyndon_min_rotation("cab") == 1);
    assert(lyndon_min_rotation("abab") == 0);
    assert(lyndon_min_rotation("baba") == 1);
    assert(lyndon_min_rotation("aaaa") == 0);
    assert(lyndon_min_rotation("dcba") == 3);
    assert(lyndon_min_rotation("z") == 0);
    assert(lyndon_min_rotation("") == 0);
    assert(lyndon_min_rotation(vector<int>{3, 1, 2, 1, 1}) == 3);

    return 0;
}
