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
    int n;
    scanf("%d", &n);
    vector<int> p(n);
    for (int& x : p) scanf("%d", &x);

    /// Plant i dies the day after every plant between it and its previous smaller one is gone
    /// best[i] is the latest death day in (prev[i], i], so hopping j = prev[j] from i - 1 covers (prev[i], i) in amortized O(1)
    auto prev = previous_smaller(p);
    vector<int> day(n, 0), best(n, 0);
    int res = 0;
    for (int i = 0; i < n; i++){
        int latest = 0;
        for (int j = i - 1; j > prev[i]; j = prev[j]) latest = max(latest, best[j]);
        day[i] = prev[i] == -1 ? 0 : latest + 1;
        best[i] = max(latest, day[i]);
        res = max(res, day[i]);
    }
    printf("%d\n", res);
    return 0;
}
