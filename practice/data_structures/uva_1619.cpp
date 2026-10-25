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
    int n, t = 0;

    while (scanf("%d", &n) == 1){
        vector<int> a(n);
        vector<long long> sum(n + 1, 0);
        for (int i = 0; i < n; i++){
            scanf("%d", &a[i]);
            sum[i + 1] = sum[i] + a[i];
        }

        auto left = previous_smaller(a), right = next_smaller(a);
        long long res = 0;
        int l = 0, r = 0;
        for (int i = 0; i < n; i++){
            int x = left[i] + 1, y = right[i] - 1;
            long long val = (sum[y + 1] - sum[x]) * a[i];
            if (val > res || (val == res && res > 0 && y - x < r - l)) res = val, l = x, r = y;
        }

        if (t++) puts("");
        printf("%lld\n%d %d\n", res, l + 1, r + 1);
    }
    return 0;
}
