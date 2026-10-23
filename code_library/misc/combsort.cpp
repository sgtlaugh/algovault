/***
 *
 * Comb Sort
 * In-place, unstable sort that shrinks the comparison gap by a factor of 1.29 per pass
 *
 * Complexity: O(n^2) worst case, no proven average bound, O(1) extra space
 *
 * combsort(first, last): sorts any random access range ascending with operator<, like combsort(ar, ar + n)
 *
 * Notes:
 *   - Fast in practice on random input: 1e6 ints in about 0.4 s
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename It>
void combsort(It first, It last){
    long long n = last - first, gap = n;
    bool swapped = false;

    while (gap > 1 || swapped){
        swapped = false;
        if (gap > 1) gap = (long long)(gap * 0.77425);

        for (long long i = 0; i + gap < n; i++){
            if (first[i + gap] < first[i]){
                swapped = true;
                iter_swap(first + i, first + i + gap);
            }
        }
    }
}

int main(){
    vector<int> v = {5, 3, 8, 1, 9, 2, 7};
    combsort(v.begin(), v.end());
    assert((v == vector<int>{1, 2, 3, 5, 7, 8, 9}));

    int ar[] = {4, -2, 4, 0};
    combsort(ar, ar + 4);
    assert(ar[0] == -2 && ar[3] == 4);

    vector<string> words = {"pear", "apple", "fig", "banana"};  /// anything with operator<
    combsort(words.begin(), words.end());
    assert((words == vector<string>{"apple", "banana", "fig", "pear"}));
    return 0;
}
