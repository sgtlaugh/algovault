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
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int m, n;
        scanf("%d %d", &m, &n);

        vector<int> h(n);
        long long res = 0;
        for (int i = 0; i < m; i++){
            int k, c;
            scanf("%d %d", &k, &c);
            for (int j = 0; k--; c ^= 1){
                int p;
                scanf("%d", &p);
                for (; p--; j++) h[j] = c ? 0 : h[j] + 1;
            }
            res = max(res, largest_rectangle(h));
        }

        printf("Case %d: %lld\n", tc, res);
    }
    return 0;
}
