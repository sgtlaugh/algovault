/***
 *
 * Main-Lorentz Tandem Repeats
 * Finds every tandem repeat (square) s[i..i + 2l) with s[i..i + l) == s[i + l..i + 2l), grouped into compact triples
 *
 * Complexity: O(n log n) time and at most n * ceil(log2 n) triples, O(n) extra memory
 *
 * tandem_repeats(s, report): calls report(first, last, l) per triple, a square s[i..i + 2l) at every first <= i <= last
 *
 * Notes:
 *   - Every tandem repeat belongs to exactly one triple, triples come in no particular order
 *   - The repeats can number Theta(n^2) ("aaaa..."): count them as the sum of last - first + 1 in a long long
 *   - Works on strings, vectors and any indexable container whose elements compare with ==, of any values
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace main_lorentz{
    /// z[i] = longest common prefix of s(i..) and s(0..), s maps an index to an element
    template <typename S>
    void z_function(S s, int n, int* z){
        if (n) z[0] = n;

        for (int i = 1, l = 0, r = 0; i < n; i++){
            int val = i < r ? min(r - i, z[i - l]) : 0;
            while (i + val < n && s(val) == s(i + val)) val++;
            if (i + val > r) l = i, r = i + val;
            z[i] = val;
        }
    }

    /// res[i] = longest common prefix of text(i..) and pattern(0..), zp = z of pattern
    template <typename P, typename S>
    void prefix_matches(P pattern, int m, const int* zp, S text, int n, int* res){
        for (int i = 0, l = 0, r = 0; i < n; i++){
            int val = i < r ? min(r - i, zp[i - l]) : 0;
            while (i + val < n && val < m && text(i + val) == pattern(val)) val++;
            if (i + val > r) l = i, r = i + val;
            res[i] = val;
        }
    }

    /// Reports the repeats inside s[lo..hi), buf holds 2 * (hi - lo) ints the children reuse as they finish first
    template <typename Container, typename F>
    void solve(const Container& s, int lo, int hi, int* buf, F& report){
        if (hi - lo < 2) return;

        int mid = lo + (hi - lo) / 2, n_left = mid - lo, n_right = hi - mid;
        solve(s, lo, mid, buf, report);
        solve(s, mid, hi, buf, report);

        /// decltype(auto) keeps const T& for ordinary containers and returns vector<bool>'s proxy by value
        auto left = [&](int i) -> decltype(auto){ return s[lo + i]; };
        auto right = [&](int i) -> decltype(auto){ return s[mid + i]; };
        auto left_reversed = [&](int i) -> decltype(auto){ return s[mid - 1 - i]; };
        auto right_reversed = [&](int i) -> decltype(auto){ return s[hi - 1 - i]; };
        int *z_left_reversed = buf, *left_vs_right = buf + n_left;
        int *z_right = buf + 2 * n_left, *right_reversed_vs_left_reversed = buf + 2 * n_left + n_right;
        z_function(left_reversed, n_left, z_left_reversed);
        z_function(right, n_right, z_right);
        prefix_matches(right, n_right, z_right, left, n_left, left_vs_right);
        prefix_matches(left_reversed, n_left, z_left_reversed, right_reversed, n_right, right_reversed_vs_left_reversed);

        /// Second half starts in left: the repeat holds the pair (mid - l, mid) and l1 >= 1 matches before mid - l
        for (int l = 1; l <= n_left; l++){
            int back = l < n_left ? z_left_reversed[l] : 0, front = left_vs_right[n_left - l];
            int lo1 = max(1, l - front), hi1 = min(l - 1, back);
            if (lo1 <= hi1) report(mid - l - hi1, mid - l - lo1, l);
        }

        /// Second half starts in right: the repeat holds pair (mid - 1, mid - 1 + l), l1 >= 1 matches up to mid - 1
        for (int l = 1; l <= n_right; l++){
            int back = right_reversed_vs_left_reversed[n_right - l], front = l < n_right ? z_right[l] : 0;
            int lo1 = max(1, l - front), hi1 = min(l, back);
            if (lo1 <= hi1) report(mid - hi1, mid - lo1, l);
        }
    }
}

template <typename Container, typename F>
void tandem_repeats(const Container& s, F report){
    int n = s.size();
    vector<int> buf(2 * n);
    main_lorentz::solve(s, 0, n, buf.data(), report);
}

template <typename F>
void tandem_repeats(const char* s, F report){
    tandem_repeats(string(s), report);
}

int main(){
    vector<pair<int, int>> found;  /// (start, half length)
    tandem_repeats("acababaee", [&](int first, int last, int l){
        for (int i = first; i <= last; i++) found.push_back({i, l});
    });
    sort(found.begin(), found.end());
    assert((found == vector<pair<int, int>>{{2, 2}, {3, 2}, {7, 1}}));  /// abab, baba and ee

    /// Count from the triples instead of listing the Theta(n^2) squares one by one
    long long total = 0;
    tandem_repeats(string(6, 'a'), [&](int first, int last, int){ total += last - first + 1; });
    assert(total == 9);  /// 5 of length 2, 3 of length 4, 1 of length 6

    int count = 0;
    tandem_repeats(vector<int>{1, 2, 1, 2, 2}, [&](int first, int last, int){ count += last - first + 1; });
    assert(count == 2);  /// 1 2 1 2 and 2 2
    return 0;
}
