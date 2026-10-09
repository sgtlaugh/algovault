/***
 *
 * Main-Lorentz Tandem Repeats
 * Finds every tandem repeat (square) s[i..i + 2l) with s[i..i + l) == s[i + l..i + 2l), grouped into compact triples
 *
 * Complexity: O(n log n) time and at most n * ceil(log2 n) triples, O(n) extra memory
 *
 * tandem_repeats(s, report) calls report(first, last, l) once per triple:
 *   s[i..i + 2l) is a tandem repeat for every first <= i <= last
 * Every tandem repeat belongs to exactly one triple, triples come in no particular order
 * The repeats themselves can number Theta(n^2) ("aaaa..."): count them as the sum of last - first + 1 in a long long
 * Works on strings, vectors and any indexable container whose elements compare with ==
 *
 * Divide and conquer: repeats crossing the middle come from z-functions of the halves and their reverses,
 * matched without separator characters, so every element value is allowed
 *
 * Example:
 *   tandem_repeats("acababaee", [&](int first, int last, int l){ ... });
 *   // repeats: abab at 2, baba at 3, ee at 7
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

    /// Reports the repeats inside s[lo..hi), those crossing mid here and the rest recursively
    /// buf holds 2 * (hi - lo) ints, reused by the children since they finish before this level writes it
    template <typename Container, typename F>
    void solve(const Container& s, int lo, int hi, int* buf, F& report){
        if (hi - lo < 2) return;

        int mid = lo + (hi - lo) / 2, nu = mid - lo, nv = hi - mid;
        solve(s, lo, mid, buf, report);
        solve(s, mid, hi, buf, report);

        /// decltype(auto) keeps const T& for ordinary containers and returns vector<bool>'s proxy by value
        auto u = [&](int i) -> decltype(auto){ return s[lo + i]; };
        auto v = [&](int i) -> decltype(auto){ return s[mid + i]; };
        auto ru = [&](int i) -> decltype(auto){ return s[mid - 1 - i]; };
        auto rv = [&](int i) -> decltype(auto){ return s[hi - 1 - i]; };
        int *zru = buf, *u_vs_v = buf + nu, *zv = buf + 2 * nu, *rv_vs_ru = buf + 2 * nu + nv;
        z_function(ru, nu, zru);
        z_function(v, nv, zv);
        prefix_matches(v, nv, zv, u, nu, u_vs_v);
        prefix_matches(ru, nu, zru, rv, nv, rv_vs_ru);

        /// Second half starts in u: the repeat holds the pair (mid - l, mid) and l1 >= 1 matches before mid - l
        for (int l = 1; l <= nu; l++){
            int back = l < nu ? zru[l] : 0, front = u_vs_v[nu - l];
            int lo1 = max(1, l - front), hi1 = min(l - 1, back);
            if (lo1 <= hi1) report(mid - l - hi1, mid - l - lo1, l);
        }

        /// Second half starts in v: the repeat holds the pair (mid - 1, mid - 1 + l) and l1 >= 1 matches up to mid - 1
        for (int l = 1; l <= nv; l++){
            int back = rv_vs_ru[nv - l], front = l < nv ? zv[l] : 0;
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
    auto squares = [](const auto& s){
        vector<pair<int, int>> res;
        tandem_repeats(s, [&](int first, int last, int l){
            for (int i = first; i <= last; i++) res.push_back({i, l});
        });
        sort(res.begin(), res.end());
        return res;
    };

    using Squares = vector<pair<int, int>>;
    assert((squares(string("")) == Squares{}));
    assert((squares(string("a")) == Squares{}));
    assert((squares(string("aa")) == Squares{{0, 1}}));
    assert((squares(string("ab")) == Squares{}));
    assert((squares(string("aaaa")) == Squares{{0, 1}, {0, 2}, {1, 1}, {2, 1}}));
    assert((squares(string("abcabczz")) == Squares{{0, 3}, {6, 1}}));
    assert((squares(string("acababaee")) == Squares{{2, 2}, {3, 2}, {7, 1}}));
    assert((squares(string("abacaba")) == Squares{}));

    assert((squares(string("##")) == Squares{{0, 1}}));
    assert((squares(string("a#a")) == Squares{}));
    assert((squares(string("a#a#")) == Squares{{0, 2}}));
    assert((squares("abaaba") == Squares{{0, 3}, {2, 1}}));
    assert((squares(vector<int>{1, 2, 1, 2, 2}) == Squares{{0, 2}, {3, 1}}));
    assert((squares(vector<bool>{1, 0, 1, 0, 0}) == Squares{{0, 2}, {3, 1}}));
    assert((squares(deque<int>{3, 1, 3, 1, 3}) == Squares{{0, 2}, {1, 2}}));

    long long total = 0;
    tandem_repeats(string(6, 'a'), [&](int first, int last, int){ total += last - first + 1; });
    assert(total == 9);

    return 0;
}
