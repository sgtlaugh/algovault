/***
 *
 * Comb Sort
 * In-place, unstable sort that shrinks the comparison gap by a factor of 1.29 per pass, https://en.wikipedia.org/wiki/Comb_sort
 *
 * Complexity: O(n^2) worst case, O(1) extra space
 * No proven average-case bound, but fast in practice on random input (1e6 ints in ~0.4 s)
 *
 * combsort(first, last) sorts any random access range ascending with operator<
 * Works on raw arrays and vectors: combsort(ar, ar + n) or combsort(v.begin(), v.end())
 * Equal elements may be reordered
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename It>
void combsort(It first, It last){
    auto n = last - first, gap = n;
    bool swapped = false;

    while (gap > 1 || swapped){
        swapped = false;
        if (gap > 1) gap = decltype(n)(gap * 0.77425);

        for (decltype(n) i = 0; i + gap < n; i++){
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

    vector<int> empty_range;
    combsort(empty_range.begin(), empty_range.end());
    assert(empty_range.empty());

    int single[] = {42};
    combsort(single, single + 1);
    assert(single[0] == 42);

    int pair_of_two[] = {2, 1};
    combsort(pair_of_two, pair_of_two + 2);
    assert(pair_of_two[0] == 1 && pair_of_two[1] == 2);

    vector<int> extremes = {INT_MAX, 0, INT_MIN, -1, INT_MAX, INT_MIN, 1};
    combsort(extremes.begin(), extremes.end());
    assert((extremes == vector<int>{INT_MIN, INT_MIN, -1, 0, 1, INT_MAX, INT_MAX}));

    vector<string> words = {"pear", "apple", "fig", "banana", "apple"};
    combsort(words.begin(), words.end());
    assert((words == vector<string>{"apple", "apple", "banana", "fig", "pear"}));

    vector<double> reals = {2.5, -1.25, 0.0, 3.75, -1.25};
    combsort(reals.begin(), reals.end());
    assert((reals == vector<double>{-1.25, -1.25, 0.0, 2.5, 3.75}));

    vector<int> reversed(1000);
    for (int i = 0; i < 1000; i++) reversed[i] = 1000 - i;
    combsort(reversed.begin(), reversed.end());
    for (int i = 0; i < 1000; i++) assert(reversed[i] == i + 1);

    mt19937 rng(20020523);
    vector<int> big(1000000);
    for (auto& x : big) x = rng() % 1000000007;
    auto expected = big;
    sort(expected.begin(), expected.end());

    clock_t start = clock();
    combsort(big.begin(), big.end());
    fprintf(stderr, "Time taken = %0.5f s\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));
    assert(big == expected);

    return 0;
}
