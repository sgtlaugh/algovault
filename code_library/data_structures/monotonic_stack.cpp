/***
 *
 * Monotonic Stack
 * Nearest smaller or greater element on each side, and the largest rectangle in a histogram
 *
 * Complexity: O(n)
 *
 * previous_smaller(a)[i]: largest j < i with a[j] < a[i], or -1
 * next_smaller(a)[i]: smallest j > i with a[j] < a[i], or n
 * previous_smaller(a, compare)[i]: largest j < i with compare(a[j], a[i]), or -1
 * next_smaller(a, compare)[i]: smallest j > i with compare(a[j], a[i]), or n
 * largest_rectangle(h): largest axis aligned rectangle under the histogram h, as long long
 *
 * compare is an ordering, strict like greater<T>() for the nearest greater element or non-strict like less_equal<T>()
 * a[i] is the minimum of every range that contains i and lies inside (previous_smaller(a)[i], next_smaller(a)[i])
 * For the sum of subarray minimums, pass less_equal<T>() to one side so a tied minimum is counted once
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
vector<int> previous_smaller(const vector<T>& a, Compare compare = Compare()){
    int n = a.size();
    vector<int> res(n), stack;
    for (int i = 0; i < n; i++){
        while (!stack.empty() && !compare(a[stack.back()], a[i])) stack.pop_back();
        res[i] = stack.empty() ? -1 : stack.back();
        stack.push_back(i);
    }
    return res;
}

template <typename T, typename Compare = less<T>>
vector<int> next_smaller(const vector<T>& a, Compare compare = Compare()){
    int n = a.size();
    vector<int> res(n), stack;
    for (int i = n - 1; i >= 0; i--){
        while (!stack.empty() && !compare(a[stack.back()], a[i])) stack.pop_back();
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
    /// Nearest strictly smaller element on each side
    vector<int> a = {5, 3, 4, 3, 1, 2, 6};
    assert((previous_smaller(a) == vector<int>{-1, -1, 1, -1, -1, 4, 5}));
    assert((next_smaller(a) == vector<int>{1, 4, 3, 4, 7, 7, 7}));

    /// Previous greater, and next smaller or equal, which stops at the tie a[3] = a[1]
    assert((previous_smaller(a, greater<int>()) == vector<int>{-1, 0, 0, 2, 3, 3, -1}));
    assert((next_smaller(a, less_equal<int>()) == vector<int>{1, 3, 3, 4, 7, 7, 7}));

    /// The bars 5 and 6 make a 2 x 5 rectangle
    assert(largest_rectangle(vector<int>{2, 1, 5, 6, 2, 3}) == 10);

    return 0;
}
