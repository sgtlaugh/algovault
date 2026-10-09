/***
 *
 * Monotonic Stack
 * Nearest strictly smaller element on each side, and the largest rectangle in a histogram
 *
 * Complexity: O(n)
 *
 * previous_smaller(a)[i]: largest j < i with a[j] < a[i], or -1
 * next_smaller(a)[i]: smallest j > i with a[j] < a[i], or n
 * So a[i] is the minimum of every range inside (previous_smaller[i], next_smaller[i]), equal values included,
 * which gives the span lengths L[i] = i - previous_smaller[i] and R[i] = next_smaller[i] - i
 * For the sum of minimums over all subarrays, sum a[i] * L[i] * R[i] counts a subarray with tied minimums
 * once per tie; make one side non-strict (stop at a[j] <= a[i]) so each subarray is credited to exactly one index
 *
 * largest_rectangle(h): largest axis aligned rectangle under the histogram h, as long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
vector<int> previous_smaller(const vector<T>& a){
    int n = a.size();
    vector<int> res(n), stack;
    for (int i = 0; i < n; i++){
        while (!stack.empty() && !(a[stack.back()] < a[i])) stack.pop_back();
        res[i] = stack.empty() ? -1 : stack.back();
        stack.push_back(i);
    }
    return res;
}

template <typename T>
vector<int> next_smaller(const vector<T>& a){
    int n = a.size();
    vector<int> res(n), stack;
    for (int i = n - 1; i >= 0; i--){
        while (!stack.empty() && !(a[stack.back()] < a[i])) stack.pop_back();
        res[i] = stack.empty() ? n : stack.back();
        stack.push_back(i);
    }
    return res;
}

template <typename T>
long long largest_rectangle(const vector<T>& h){
    auto left = previous_smaller(h), right = next_smaller(h);
    long long best = 0;
    for (int i = 0; i < (int)h.size(); i++) best = max(best, (long long)h[i] * (right[i] - left[i] - 1));
    return best;
}

int main(){
    vector<int> a = {5, 3, 4, 3, 1, 2, 6};
    auto left = previous_smaller(a), right = next_smaller(a);
    assert((left == vector<int>{-1, -1, 1, -1, -1, 4, 5}));
    assert((right == vector<int>{1, 4, 3, 4, 7, 7, 7}));

    vector<int> L(7), R(7);
    for (int i = 0; i < 7; i++) L[i] = i - left[i], R[i] = right[i] - i;
    assert((L == vector<int>{1, 2, 1, 4, 5, 1, 1}));
    assert((R == vector<int>{1, 3, 1, 1, 3, 2, 1}));

    assert(largest_rectangle(vector<int>{2, 1, 5, 6, 2, 3}) == 10);
    assert(largest_rectangle(vector<int>{}) == 0);
    assert(largest_rectangle(vector<int>{4}) == 4);
    assert(largest_rectangle(vector<long long>{1000000000, 1000000000, 1000000000}) == 3000000000LL);
    return 0;
}
