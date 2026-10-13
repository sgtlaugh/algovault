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
    int p;
    scanf("%d", &p);
    while (p--){
        int s, b;
        scanf("%d %d", &s, &b);
        vector<vector<bool>> occupied(s, vector<bool>(s, false));
        while (b--){
            int r1, c1, r2, c2;
            scanf("%d %d %d %d", &r1, &c1, &r2, &c2);
            for (int i = r1 - 1; i < r2; i++){
                for (int j = c1 - 1; j < c2; j++) occupied[i][j] = true;
            }
        }

        /// Row by row, the free run above each square is a histogram whose largest rectangle ends on that row
        vector<int> height(s, 0);
        long long res = 0;
        for (int i = 0; i < s; i++){
            for (int j = 0; j < s; j++) height[j] = occupied[i][j] ? 0 : height[j] + 1;
            res = max(res, largest_rectangle(height));
        }
        printf("%lld\n", res);
    }
    return 0;
}
